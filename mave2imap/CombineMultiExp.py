#!/usr/bin/env python3

import lr_tools as lt

import sys
import os
import statistics
import re

from datetime import datetime
import time
import tempfile
import subprocess
import argparse
import pickle
import pprint
import fileinput
from copy import deepcopy
from mpi4py import MPI
from itertools import islice

import fileinput as FI
import h5py

try:
    import configparser
except:
    import ConfigParser as configparser
from configparser import ExtendedInterpolation



class Process:
    def __init__(self, rootname_before, rootname_after, working_dir, count_threshold, nameout, mpi, config_file, pdbfile=False, verbose = False, cleanup = False):
        self.nbef = rootname_before
        self.naft = rootname_after
        self.working_dir = working_dir

        self.nameout = nameout
        self.mpi = mpi
        self.config_file = config_file
        self.log = ""
        self.pdbfile = pdbfile
        self.verbose = verbose
        self.cleanup = cleanup

        # Parameters to pass as options
        self.use_MONO_for_statistics_onWTsilMut = True
        self.use_MONO_for_statistics_onSynMut = True
        self.use_distinct_stats_between_MONO_vs_DITRI = True

        # Analyze count threshold and create a list of threshold for every categories of mutant, MONO, DI or TRI
        if count_threshold:
            if count_threshold == "sum_all":
                self.count_threshold = "sum_all"
            else:
                nthresh = count_threshold.split(',')
                if len(nthresh) == 1:
                    self.count_threshold = [int(count_threshold),int(count_threshold),int(count_threshold)]
                else:
                    self.count_threshold = [int(i) for i in nthresh]
        else:
            self.count_threshold = None


        self.read_config_file()

        self.d3to1, self.d1to3 = lt.get_dicaa()
        self.dmut = lt.parse_twist_devis_table(self.mutfile, self.dn2a, self.d3to1, start=None, stop=None,
                                               info_codon_index=None)
        self.d_enrich = dict()
        self.run_process()

        if self.pdbfile:
            self.map_on_pdb()


    def run_process(self):

        os.chdir(self.working_dir)
        self.comparecond_out = "{}_compare_conditions.out".format(self.nameout)
        self.dmut_cond = {}
        self.mutant_out = "{}_mutations.out".format(self.nameout)
        self.dbcode = {}
        self.dmutcode = {}
        self.dmutcode['order'] = []

        # Dic of all WT silence nuc mutations with counts >0
        self.dwtsilent = {}
        # Dic of all pairs of synonymous mutations with counts >0
        self.dsyno = {}

        for fname_condition in [self.nbef, self.naft]:
            print("Analyzing condition ",fname_condition)
            fbcode = '{}_assembled.out'.format(fname_condition)
            path_fcode = os.path.join(self.working_dir,fbcode)
            mutcode = '{}_mutants_for_rev.out'.format(fname_condition)
            path_mutcode = os.path.join(self.working_dir,mutcode)

            # Analyze the mutant counts for every primer output
            with (open(path_mutcode)) as fin:
                for line in fin:
                    s = line.split()
                    nucmut = s[0]
                    mutation_type = len(nucmut.split('.')) - 1 # 0, 1 or 2 for MONO, DI or TRI
                    aamut = s[1]
                    aalist = [int(s[1][1:-1]),s[1][-1]] # index de mutation, AA mutant
                    if not self.count_threshold:
                        count = eval(s[2])
                        distrib = s[3]
                    elif self.count_threshold == "sum_all":
                        # We sum the entire number of rerads with a given mutation, irrespective of the barcode
                        distrib = s[3]
                        count = sum([int(kk) for kk in distrib.split('_')])
                    else:
                        count, distrib = self.recompute_count_with_thresh(s[3], self.count_threshold[mutation_type])
                    if nucmut not in self.dmutcode:
                        self.dmutcode[nucmut] = {}
                        self.dmutcode['order'].append([aalist[0],aalist[1],nucmut]) # index de mutation AA, AA mutant, nucmut
                    self.dmutcode[nucmut][fname_condition] = [aamut,count,distrib]
                    #print(self.dmutcode[nucmut][fname_condition])
        self.dmutcode['order'].sort()
        #print(self.dmutcode['order'])

        # Get the count values from WT silents
        numberOf_WTenrich, median_WTenrich = self.analyze_silentWT(self.nbef,self.naft)
        print(median_WTenrich)
        # Compare the enrichments between synonymous mutations
        numberOf_Syno, syno_Mean, syno_StdDev = self.analyze_synonymous(self.nbef, self.naft, median_WTenrich)

        self.normalization_factor = median_WTenrich
        self.perturbation_threshold = syno_StdDev

        # Dump the counts into a file
        with(open(self.comparecond_out,"w")) as fout:
            fout.write(self.log)
            fout.write("# COLUMNS = 1) Nucleotide Mutant\n"
                       "#           2) AminoAcid Mutant\n"
                       "#           3) Mutant counts before selection \n"
                       "#           4) Mutant counts after selection \n"
                       "#           5) Enrichment corrected using WTsilent and Synonyms \n")
            fout.write("#\n")

            for mutlist in self.dmutcode['order']:
                mut = mutlist[-1]
                number_of_mutations = len(nucmut.split('.'))
                if number_of_mutations == 1:
                    mutant_category = 0
                else:
                    mutant_category = 1
                if self.use_distinct_stats_between_MONO_vs_DITRI:
                    try:
                        normalization_factor = self.normalization_factor[mutant_category]
                    except IndexError:
                        normalization_factor = self.normalization_factor[0]
                aamut = self.dmutcode[mut][self.nbef][0]
                fout.write("{}\t{}\t".format(mut, aamut))
                dcount = {}
                for fname_condition in [self.nbef, self.naft]:
                    dcount[fname_condition] = self.dmutcode[mut][fname_condition][1]
                    fout.write("{}\t".format(dcount[fname_condition]))
                if dcount[self.nbef] != 0:
                    enrich = dcount[self.naft] * (1./normalization_factor) / dcount[self.nbef]
                    self.d_enrich[mut] = enrich
                    fout.write("{:.2f}\n".format(enrich))
                else:
                    enrich = "Nan"
                    fout.write("{}\n".format(enrich))

    def map_on_pdb(self):
        d_aa2enrich = dict()
        d_aa2pert  = dict()
        d_aa2enrich['order'] = []
        self.pdb_out = "{}_map_perturbation.pdb".format(self.nameout)
        for l_mut in self.dmutcode['order']: # index de mutation AA, AA mutant, nucmut
            if l_mut[-1] not in self.d_enrich:
                continue
            if l_mut[-1] not in self.dmut:
                continue
            enrich = self.d_enrich[l_mut[-1]]
            aamut_index = "{}{}".format(l_mut[1],l_mut[0]) # example A50
            if not aamut_index in d_aa2enrich:
                d_aa2enrich[aamut_index] = []
                d_aa2enrich['order'].append(aamut_index)

            d_aa2enrich[aamut_index].append(enrich)
            print(l_mut,enrich)
        print("perturbation_threshold:",self.perturbation_threshold)
        for aamut in d_aa2enrich['order']:
            mean_enrich =  sum(d_aa2enrich[aamut])*1./len(d_aa2enrich[aamut])
            low_thresh = 1.- self.perturbation_threshold
            high_thresh = 1.+ self.perturbation_threshold
            print(aamut,mean_enrich)
            if aamut[0]=="Z":
                continue
            if mean_enrich > low_thresh and mean_enrich < high_thresh:
                perturbation = 0.
            else:
                perturbation = min(abs(low_thresh-mean_enrich),abs(high_thresh-mean_enrich))
            index = aamut[1:]
            """
            if index=='65' or index=='71': # Excetion for VHH1
                continue
            """
            if perturbation == 0.:
                continue
            if index not in d_aa2pert:
                d_aa2pert[index] = 0
            d_aa2pert[index] += perturbation

        lt.map_perturbation_score2pdb(self.pdbfile, d_aa2pert, self.pdb_out)#,chain_pdb=["G"])

    def analyze_silentWT(self,before,after):
        """
        Return the median of the enrichments for WT silent mutations
        :param before: root_name of the file with mutant counts before selection
        :param after: root_name of the file with mutant counts after selection
        :return:
        """
        # Get the count values from WT silents
        for index,mutaa,nucmut in self.dmutcode['order']:
            number_of_mutations = len(nucmut.split('.'))
            if number_of_mutations == 1:
                mutant_category = 0
            else:
                mutant_category = 1
            if not self.use_MONO_for_statistics_onWTsilMut and number_of_mutations == 1: # We don't consider single mutations due to potential noise
                continue
            aamut,count_bef,distrib = self.dmutcode[nucmut][before]
            aamut,count_aft,distrib = self.dmutcode[nucmut][after]
            #print(aamut, mutaa)
            if aamut[0] == aamut[-1] and count_bef > 0:
                enrich = count_aft*1. / count_bef
                self.dwtsilent[nucmut] = [count_bef, count_aft, enrich, mutant_category]
                #print("WT: {}, {}, {}, {:.3f}".format(nucmut, count_bef, count_aft, enrich))

        list_cases = ["MONO", "DI+TRI"]
        if self.use_distinct_stats_between_MONO_vs_DITRI:
            list_WTenrich = [[] for ii in range(2)]
            list_WTenrich[0] = [self.dwtsilent[i][2] for i in self.dwtsilent if self.dwtsilent[i][3]==0]
            list_WTenrich[1] = [self.dwtsilent[i][2] for i in self.dwtsilent if self.dwtsilent[i][3]==1]
        else:
            list_WTenrich = [[self.dwtsilent[i][2] for i in self.dwtsilent]]
        median_WTenrich = []
        numberOf_WTenrich = []
        cases = []
        for jj in range(len(list_WTenrich)):

            if len(list_WTenrich[jj]) == 0:
                continue
            median_WTenrich.append(statistics.median(list_WTenrich[jj]))
            numberOf_WTenrich.append(len(list_WTenrich[jj]))
            cases.append(list_cases[jj])

        if len(numberOf_WTenrich) > 1:
            message = "{} silent WT mutants {} have a median ratio of counts before and after at {} {}\n"\
                .format(numberOf_WTenrich, cases, median_WTenrich, cases)
        else:
            message = "{} silent WT mutants {} have a median ratio of counts before and after  at {}\n"\
                .format(numberOf_WTenrich, cases, median_WTenrich)
        print("> {}".format(message))
        self.log += "# {}".format(message)

        message += "Ratio of counts will be used as normalization factor to center enrichment values around 1\n"
        print("> {}".format(message))
        self.log += "# {}".format(message)

        return numberOf_WTenrich,median_WTenrich

    def analyze_synonymous(self,before, after, normalization_factor):
        """
        Analyze the standard deviation between the enrichments calculated between synonymous pairs.
        :param before: root_name of the file with mutant counts before selection
        :param after: root_name of the file with mutant counts after selection
        :param normalization: Normalization factor used to compare after and before selections
        :return: standard deviation of the differences between the synonymous pairs.
        """
        # oscar
        print(before, after, normalization_factor)
        set_nucmut_analyzed = set() # set of mutants already seen was soi as not to repeat analysis in symetry

        for index,mutaa,nucmut in self.dmutcode['order']:
            number_of_mutations = len(nucmut.split('.'))
            if number_of_mutations == 1:
                mutant_category = 0
            else:
                mutant_category = 1
            if self.use_distinct_stats_between_MONO_vs_DITRI:
                try:
                    float_normalization_factor = normalization_factor[mutant_category]
                except IndexError:
                    float_normalization_factor = normalization_factor[0]
            else:
                float_normalization_factor = normalization_factor[0]
            if not self.use_MONO_for_statistics_onSynMut and number_of_mutations == 1: # We don't consider single mutations due to potential noise
                continue
            aamut,count_bef,distrib = self.dmutcode[nucmut][before]
            aamut,count_aft,distrib = self.dmutcode[nucmut][after]
            if not self.dmut[nucmut]['nuc_syn']:
                continue
            if count_bef == 0:
                continue
            enrich = count_aft * (1./float_normalization_factor) / count_bef
            for syno in self.dmut[nucmut]['nuc_syn']:
                if not self.use_MONO_for_statistics_onSynMut and number_of_mutations == 1: # We don't consider single mutations due to potential noise
                    continue
                if syno in set_nucmut_analyzed:
                    continue
                number_of_mut_in_syno = len(syno.split('.'))
                if number_of_mut_in_syno == 1:
                    syno_category = 0
                else:
                    syno_category = 1
                if self.use_distinct_stats_between_MONO_vs_DITRI:
                    try:
                        float_normalization_factor2 = normalization_factor[mutant_category]
                    except IndexError:
                        float_normalization_factor2 = normalization_factor[0]
                else:
                    float_normalization_factor2 = normalization_factor[0]
                aamut_syn,count_bef_syn,distrib_syn = self.dmutcode[syno][before]
                aamut_syn,count_aft_syn,distrib_syn = self.dmutcode[syno][after]
                if count_bef_syn == 0:
                    continue
                enrich_syno = count_aft_syn * (1./float_normalization_factor2) / count_bef_syn
                pairsyno_name = "{}_{}&{}".format(aamut,nucmut,syno)
                diff_enrich = enrich - enrich_syno
                self.dsyno[pairsyno_name] = [enrich,enrich_syno,diff_enrich]
            set_nucmut_analyzed.add(nucmut)
        list_diffEnrich = [self.dsyno[i][2] for i in self.dsyno]
        list_couples = [(self.dsyno[i][0],self.dsyno[i][1]) for i in self.dsyno]
        syno_StdDev = statistics.stdev(list_diffEnrich)
        syno_Mean = statistics.mean(list_diffEnrich)
        syno_Length = len(list_diffEnrich)

        message = "{} pairs of synonymous mutations were analyzed\n".format(syno_Length)
        print("> {}".format(message))
        self.log += "# {}".format(message)

        message = "Mean value of the difference between enrichments is : {:2f}\n".format(syno_Mean)
        print("> {}".format(message))
        self.log += "# {}".format(message)

        message = "StdDev value of the difference between enrichments is : {:2f}\n".format(syno_StdDev)
        print("> {}".format(message))
        self.log += "# {}".format(message)

        #pprint.pprint(list_couples)

        return syno_Length,syno_Mean,syno_StdDev

    def recompute_count_with_thresh(self, distrib, thresh):
        list_size = [int(size) for size in distrib.split('_')]
        new_distrib = []
        for elt in list_size:
            if elt >= thresh:
                new_distrib.append(elt)
        return len(new_distrib), '_'.join([str(i) for i in new_distrib])

    def read_config_file(self):
        """

        :return:
        """
        config = configparser.ConfigParser(interpolation=ExtendedInterpolation())
        config.read(self.config_file)

        self.exe_dir    = config.get("Paths", "EXE_DIR")
        self.qqsub_exe  = config.get("Paths", "QSUB_EXE")

        self.use_qqsub   = eval(config.get("Parameters", "USE_QSUB"))

        self.mutfile       = config.get("Files", "TWIST_MUT_TABLE")
        faa2nuc            = config.get("Files", "AA2NUC")
        self.da2n = pickle.load(open(faa2nuc,'rb'),encoding='bytes') # due to py2 to 3
        self.da2n = lt.convert_dic_bytes2string(self.da2n) # due to py2 to 3
        fnuc2aa   = config.get("Files", "NUC2AA")
        self.dn2a = pickle.load(open(fnuc2aa,'rb'))#,encoding='bytes')


def Main():

    try:
        default_configfile = os.path.join(os.path.dirname(os.path.abspath(os.readlink(os.path.abspath(__file__)))),
                                  "NGS_DMS.ini")
    except:
        default_configfile = os.path.join(os.path.dirname(os.path.abspath(os.path.realpath(os.path.abspath(__file__)))),
                                  "NGS_DMS.ini")


    usage = "Script to parse the reads from a long read experiment.\n" + \
            "output is a dictionary encoded as a hdf5 file"
    parser = argparse.ArgumentParser(usage)
    parser.add_argument('-b', '--before_selection', action="store", help='root name <>_assembled.out + <>_mutants_for_rev.out before selection ')
    parser.add_argument('-a', '--after_selection', action="store", help='root name <>_assembled.out + <>_mutants_for_rev.out after selection')
    parser.add_argument('-p', '--pdb_reference', action="store", default=False, help='pdbfile used to map the perturbations')
    parser.add_argument('-d', '--working_directory', action="store", default='./', help='Working directory')
    parser.add_argument('-t', '--threshold', action="store", default=None, help='Which threshold to use for minimal number of counts. '
                                                                                'Default just takes the counts how they are provided')
    parser.add_argument('-o', '--nameout', action="store", default='output', help='output_file')
    parser.add_argument('--mpi', action="store_true", default=False, help='triggers mpi running')
    parser.add_argument('-c', '--config_file', default=default_configfile, action="store", help='name of the configuration file')
    parser.add_argument('--verbose', action="store_true", default=False, help="print extra info")
    parser.add_argument('--cleanup', action="store_true", default=False, help="cleanup temporary files after execution")

    if len(sys.argv) == 1:
        print("type -h or --help for more help information")
        sys.exit(1)

    args = parser.parse_args()

    if args.config_file == default_configfile:
        os.system("cp {} {}".format(default_configfile,os.getcwd()))

    if args.before_selection and args.after_selection :
        """
        Here generate 
        """
        P = Process(args.before_selection,
                      args.after_selection,
                      args.working_directory,
                      args.threshold,
                      args.nameout,
                      args.mpi,
                      args.config_file,
                      pdbfile = args.pdb_reference,
                      verbose = args.verbose,
                      cleanup = args.cleanup)

    else:
        print("type -h or --help for more help information")
        sys.exit(1)

if __name__ == "__main__":


    Main()
