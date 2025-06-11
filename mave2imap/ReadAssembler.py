#!/usr/bin/env python3

import lr_tools as lt

import sys
import os
import re

from datetime import datetime
import time
import tempfile
import argparse
import pickle
import fileinput
from pprint import pprint
import copy
from mpi4py import MPI
from itertools import islice

import fileinput as FI
import h5py

try:
    import configparser
except:
    import ConfigParser as configparser
from configparser import ExtendedInterpolation



class Assembler:
    def __init__(self, fforward , freverse, fileout, config_file, quality_report, log_file, tag,
                 use_mpi= False, tandem_barcode = False, verbose = False, cleanup = True):
        self.fforward = fforward
        self.freverse = freverse
        self.root_fname = self.get_root_filename(fileout)
        self.config_file = config_file
        self.verbose  = verbose
        self.cleanup  = cleanup
        self.tandem_barcode = tandem_barcode
        self.use_mpi = use_mpi
        if tag:
            self.tag = "_" + tag
        else:
            self.tag = ""

        self.dbarcode = None

        self.thresh_rejection_number_of_observations_per_barcode = 0 # Keep this to 0 to have a fair estimate of
                                                                # the results in case no barcode would have been used.
        self.thresh_maxnumber_of_variants_per_barcode = 6
        self.thresh_number_of_observations_in_multimut_per_barcode = 2
        self.max_interbc_dist = 14 # max distance separating 2 barcodes (7+7, X does not count, dephaser diff => max_dist
        self.thresh_abundance = [(1, 2), (2, 3), (3, 50), (50, 10000)]  # [(1,2)(2,3)(4,10)(10,50)(50,10000)]

        self.quality_report = quality_report
        if self.quality_report:
            print(">>> Quality report activated. Quality report files will be generated.")
            self.dquality = {}
            self.dqcode = {}
            for category in ["W","T","D","S","P"]:
                self.dquality[category] = {}
                self.dquality[category]["low"] = []
                self.dquality[category]["medium"] = []
                self.dquality[category]["high"] = []
                self.dqcode[1] = "S"
                self.dqcode[2] = "D"
                self.dqcode[3] = "T"

        self.d3to1, self.d1to3 = lt.get_dicaa()

        self.start_mpi()

        if self.rank == 0:
            print("\n>>> Reading the configuration file.")

        self.read_config_file()
        self.dmut = self.get_dic_expect_mutations()

        self.run_process()

    def run_process(self):

        #
        # Each read hexad code is assigned a key as a list of properties for each read
        #   [hexa_code] = [variant_code, quality, in_mut_table, barcode]
        #
        number_of_lines_freverse = self.get_nlines(self.freverse)
        number_of_lines_fforward = self.get_nlines(self.fforward)

        self.dreverse,self.lreverse = self.parse_file(self.freverse, number_of_lines_freverse, get_barcode = True)
        if self.rank == 0:
            print("> Size of reads dreverse on every working cpu:",len(self.dreverse))

        self.dforward,self.lforward = self.parse_file(self.fforward, number_of_lines_fforward, get_barcode = self.tandem_barcode)
        if self.rank == 0:
            print("> Size of reads dforward on every working cpu:",len(self.dforward))

        """
        # For NGS_DMS, we need to group together the barcode and decipher the nature of the mutation
        """
        if self.rank==0:
            print(">>> Running self.match_mutant_to_barcodes()")
        dbarcode = self.match_mutant_to_barcodes(fileout = "distrib_barcodes_vs_global_assemby_cases{}.out".format(self.tag))
                # Content of dbarcode key:barcode value:dict(mutant:number_of_reads)
        if self.rank ==0:
            print("> Size of dic of dbarcode:",len(dbarcode))

        # flush memory
        self.dforward = None
        self.dreverse = None
        self.lreverse = None
        self.lforward = None

        dbc_cat = self.sort_barcodes_categories(dbarcode, fileout = "distrib_barcodes_vs_mutant_types{}.out".format(self.tag))

        self.print_to_file_distrib_bc_interdistance(dbc_cat,dir_distrib="distrib_barcode_similarities{}".format(self.tag), fileout="distrib_bcdist{}".format(self.tag))

        # The end of the script runs without mpi

        if self.rank == 0:
            self.print_to_file_compare_bad_correct_mono(dbc_cat, dbarcode, "distrib_bad_mono_in_seq{}.out".format(self.tag), length_protein=(len(self.coding_seq)//3)+2)

            # Eventually not found useful to improve the signal-to-noise ratio
            if self.quality_report:
                print(self.dquality)
                self.analyze_quality()

            self.dobsmutants = self.count_mutants(self.root_fname,dbarcode,code_forw_rev='assembled')
            print("> Number of mutants analyzed :",len(self.dobsmutants))

            # Dump the mutant counts in a file mutant_rev_for.out
            self.analyze_master_mutant()

            """#################
             END OF MAIN SCRIPT
            """#################

    def start_mpi(self):
        if self.use_mpi:
            # MPI communication
            self.comm = MPI.COMM_WORLD
            self.rank = self.comm.Get_rank()
            self.size = self.comm.Get_size()
        else:
            self.comm = None
            self.rank = 0
            self.size = 1

    def get_nlines(self,file):
        if self.rank == 0:
            number_of_lines = lt.count_file_lines(file)
            #print("I finished")
        else:
            number_of_lines = None
        if self.use_mpi:
            self.comm.Barrier()
            number_of_lines = self.comm.bcast(number_of_lines, root=0)
        if self.rank==0:
            print("> The number of reads in the input file is: {}\n".format(number_of_lines))
        return number_of_lines

    def get_root_filename(self,fileout):
        """

        :return:
        """
        s = fileout.split('.')
        if len(s) == 1:
            return fileout
        else:
            return '.'.join(fileout.split('.')[:-1])

    def get_dic_expect_mutations(self):
        """

        :return:
        """
        return lt.parse_twist_devis_table(self.mutfile, self.dn2a, self.d3to1, start=None, stop=None, info_codon_index=None)

    def parse_file(self,file, number_of_lines, get_barcode = False):
        """
        File parsing. Can be used using mpi. Every cpu is taking charge of a chunk of the file.
        This chunk will be passed to the next method so that cpu will work independently on its own chunk
        Only after forward and reverse sets of reads are fused, the barcodes are extracted and mutant assigned.
        The program will then work on the whole library of reads.

        :param file: a one-line for every read file starting with read hexadecimal code, the quality score,
        the nucleotide mutant code and the amino acid mutant code:
        "0x2540be402	0.13	N	12A_42A.43T.44T_46A_100A	AACACAGC"
        Forward and reverse file should  have exactely the same number of lines.
        Since indexing of reads in the file should match exactely between the forward and the reverse file.

        :param get_barcode: extracts the barcode when available.

        :return: Every cpu returns a list of read codes and a dictionary with the following info:
        d[read code] = [variant_code, quality, in_mut_table, barcode]
        """
        # Define the start and end of file reading for every cpu
        size_chunks = (number_of_lines) / self.size + 1 #
        index_start = int(self.rank * (size_chunks))
        if self.rank == self.size - 1:
            index_stop = (number_of_lines)
        else:
            index_stop = int((self.rank + 1) * (size_chunks))

        # Initialization of list and dictionary for every cpu
        if self.use_mpi:
            self.comm.Barrier()
            list_dread2bc = [dict() for i in range(self.size)]
            list_lread2bc = [list() for i in range(self.size)]
            d_read2barcodemut = self.comm.scatter(list_dread2bc,root=0)
            l_read2barcodemut = self.comm.scatter(list_lread2bc,root=0)
            if self.rank==0:
                print(">>> Parsing the file {} using {} cpus".format(file, self.size))
        else:
            print(">>> Parsing the file {} ".format(file))
            d_read2barcodemut = dict()
            l_read2barcodemut = list()

        with open(file) as file_obj:
            n=1
            object_to_iterate = islice(file_obj, index_start, index_stop)
            for line in object_to_iterate:
                if line[0] == "#":
                    continue
                s = line.split()
                read_hexidx  = int(s[0],16)
                quality = int(eval(s[1])*100)
                variant_code = s[3]
                in_mut_table = s[2]
                if get_barcode:
                    try:
                        barcode = s[4]
                    except IndexError:
                        barcode = None
                    d_read2barcodemut[read_hexidx] = [variant_code, quality, in_mut_table, barcode]
                    l_read2barcodemut.append(read_hexidx)
                else:
                    d_read2barcodemut[read_hexidx] = [variant_code, quality, in_mut_table]
                    l_read2barcodemut.append(read_hexidx)
        if self.use_mpi:
            self.comm.Barrier()
            # We don't gather the dictionaries. We return every chunk independently
            """
            d_read2barcodemut = self.comm.gather(d_read2barcodemut, root=0)
            l_read2barcodemut = self.comm.gather(l_read2barcodemut, root=0)
            if self.rank == 0:
                print(">>> Gathered back a list of {} dictionary.".format(self.size))
                for dd in d_read2barcodemut[1:]:
                    d_read2barcodemut[0] = lt.data_merge(d_read2barcodemut[0],dd)
                d_read2barcodemut = d_read2barcodemut[0]
                for ll in l_read2barcodemut[1:]:
                    l_read2barcodemut[0] = lt.data_merge(l_read2barcodemut[0],ll)
                l_read2barcodemut = l_read2barcodemut[0]
                print(">>> Merged them into a single dictionary")
            else:
                l_read2barcodemut = None
                d_read2barcodemut = None
            l_read2barcodemut = self.comm.bcast(l_read2barcodemut, root=0)
            d_read2barcodemut = self.comm.bcast(d_read2barcodemut, root=0)
            if self.rank==0:
                print(">>> Result Broadcast to all cpus")
            """
        return d_read2barcodemut, l_read2barcodemut

    def analyze_quality(self):
        """
        Test to integrate qualities to filter reads but was not kept as an interesting method eventually
        :param tandem_barcodes:
        :return:
        """
        for k in self.dreverse:
            if k in self.dforward:
                forward = self.dforward[k]
                reverse = self.dreverse[k]
                barcode_fow = forward[-1]
                barcode_rev = reverse[-1]
                quality_fow = int(eval(forward[1]) * 100)
                quality_rev = int(eval(reverse[1]) * 100)
                if not (barcode_fow and barcode_rev):  # Both forward and reverse are None
                    self.dquality["P"]["medium"].append(quality_fow)
                    self.dquality["P"]["medium"].append(quality_rev)
                    continue
                # Foward and Reverse Barcodes are merged with an X in the middle
                full_barcode = barcode_fow + "X" + barcode_rev

                mcode_fow = forward[2]
                mcode_rev = reverse[2]

                mutant_fow = forward[0]
                mutant_rev = reverse[0]
                if not full_barcode in self.dquality:
                    self.dquality[full_barcode] = {}
                if mcode_fow + mcode_rev == "OO":
                    # Case 1 : Full agreement between both reads on an expected mutation
                    if mutant_fow == mutant_rev:
                        if not mutant_fow in self.dquality[full_barcode]:
                            self.dquality[full_barcode][mutant_fow] = []
                        self.dquality[full_barcode][mutant_fow].append(quality_fow)
                        self.dquality[full_barcode][mutant_fow].append(quality_rev)
                elif mcode_fow + mcode_rev == "NN":
                    # Case 4: both reads are WT
                    if mutant_fow + mutant_rev == "WW":
                        if not "W" in self.dquality[full_barcode]:
                            self.dquality[full_barcode]["W"] = []
                        self.dquality[full_barcode]["W"].append(quality_fow)
                        self.dquality[full_barcode]["W"].append(quality_rev)
                if len(self.dquality[full_barcode]) == 0:
                    del self.dquality[full_barcode]
            else:
                if self.verbose:
                    print("READ {} NOT FOUND in forward file".format(k))

    def match_mutant_to_barcodes(self, fileout=None, tandem_barcodes = False):
        """
        This method analyzes all the hexacodes for every lines, fuses the barcodes
        and decides the status of the mutation
        Distinguishes 12 categories of mutations

        :param tandem_barcodes: used in a previous version of the script for long read analysis
        :return: a dictionary of barcodes with information about the mutants associated with it.
        """

        if self.rank==0:
            print(">>> Matching forward and reverse reads. Counting barcodes and associated mutations")

        self.list_attributes = ['Expect_mut_agree','WT_agree','WT+Expect_mut',
                                'Expect_DITRI+Bullshit','Expect_mut_disagree_DITRIwin', 'Expect_MONOmut_disagree',
                                'Expect_DITRImut_disagree','ConstantTail_issue','NotExpected',
                                'NotExpect_ContainDIorTRI_retrieved','NotExpect_ContainDIorTRI_notkept',
                                'NotExpect_MONO_kept','NotExpect_DI_kept','NotExpect_TRI_kept']
        dtranslate = {key:value for key, value in [(1,"MONO"),(2,"DI"),(3,"TRI")]}

        len_freads = len(self.lforward)
        len_rreads = len(self.lreverse)
        if len_freads != len_rreads:
            print("ERROR: the number of reads in the forward and the reverse set of reads are not equal. "
                  "Process is stopped."
                  "Check what's going wrong first.")
            sys.exit()

        # Define the start and end of read list analysis for every cpu
        size_chunks = (len_rreads) / self.size + 1 #
        index_start = int(self.rank * (size_chunks))
        if self.rank == self.size - 1:
            index_stop = len_rreads
        else:
            index_stop = int((self.rank + 1) * (size_chunks))

        # Initialisation of the dict
        ### MPI block
        if self.use_mpi:
            self.comm.Barrier()
            if self.rank == 0:
                dbarcode = [dict() for i in range(self.size)]
                count_cases = [{key:0 for key in self.list_attributes} for i in range(self.size)] # Used to perform some statistics
                number_of_common_hexindexes = [0 for i in range(self.size)]
                number_of_unmatched_hexindexes = [0 for i in range(self.size)]
            else:
                dbarcode = None
                count_cases = None
                number_of_common_hexindexes = None
                number_of_unmatched_hexindexes = None
            dbarcode = self.comm.scatter(dbarcode, root=0)
            count_cases = self.comm.scatter(count_cases, root=0)
            number_of_common_hexindexes = self.comm.scatter( number_of_common_hexindexes, root=0)
            number_of_unmatched_hexindexes = self.comm.scatter( number_of_unmatched_hexindexes, root=0)
            if self.rank==0:
                print(">>> Scattered a list of {} dictionary of barcodes and count_cases in different mutant categories".format(self.size))
        else:
            dbarcode = dict()
            count_cases = {key:0 for key in self.list_attributes}
            number_of_common_hexindexes = 0
            number_of_unmatched_hexindexes = 0

        for ii, reverse_code in enumerate(self.lreverse):
            forward_code = self.lforward[ii]
            forward = self.dforward[forward_code] # a list of elements of type d[read code] = [variant_code, quality, in_mut_table, barcode]
            reverse = self.dreverse[reverse_code]
            number_of_common_hexindexes += 1
            barcode_fow = forward[-1]
            barcode_rev = reverse[-1]
            if not (barcode_fow and barcode_rev):  # Either forward or reverse are None
                count_cases['ConstantTail_issue'] += 1
                continue
            # Foward and Reverse Barcodes are merged with an X in the middle
            full_barcode = barcode_fow + "X" + barcode_rev

            # Code is to know if the mutation is canonical or not. Can be O or N
            mcode_fow = forward[2]
            mcode_rev = reverse[2]
            # mutant is the nucleotide mutation
            mutant_fow = forward[0]
            mutant_rev = reverse[0]

            if not full_barcode in dbarcode:
                dbarcode[full_barcode] = {}
            if mcode_fow + mcode_rev == "OO":
                # Case 1 : Full agreement between both reads on an expected mutation
                if mutant_fow == mutant_rev:
                    if not mutant_fow in dbarcode[full_barcode]:
                        dbarcode[full_barcode][mutant_fow] = 0
                    dbarcode[full_barcode][mutant_fow] += 1
                    count_cases['Expect_mut_agree'] += 1
                # Case 2 : Disagreement between expected mutation...
                #          We keep it in the report but not in the enrichments
                else:
                    # Case 2a: One mutation is a di- or tri_ nucleotide variant. The other a mono
                    # Our decision = select the longest variant
                    if len(mutant_fow.split('.')) >= 2 and len(mutant_rev.split('.')) == 1:
                        if not mutant_fow in dbarcode[full_barcode]:
                            dbarcode[full_barcode][mutant_fow] = 0
                        dbarcode[full_barcode][mutant_fow] += 1
                        count_cases['Expect_mut_disagree_DITRIwin'] += 1
                    # Case 2b: One mutation is a di- or tri_ nucleotide variant. The other a mono
                    elif len(mutant_rev.split('.')) >= 2 and len(mutant_fow.split('.')) == 1:
                        if not mutant_rev in dbarcode[full_barcode]:
                            dbarcode[full_barcode][mutant_rev] = 0
                        dbarcode[full_barcode][mutant_rev] += 1
                        count_cases['Expect_mut_disagree_DITRIwin'] += 1
                    # Case 2c: Both mutations are di- or tri_ nucleotide variant. We don't keep.
                    elif len(mutant_rev.split('.')) >= 2 and len(mutant_fow.split('.')) >= 2:
                        mutation = "_".join([mutant_fow, mutant_rev])
                        if not mutation in dbarcode[full_barcode]:
                            dbarcode[full_barcode][mutation] = 0
                        dbarcode[full_barcode][mutation] += 1
                        count_cases['Expect_DITRImut_disagree'] += 1
                    # Case 2d:  Both mutations are mono nucleotide variant. We don't keep.
                    else:
                        mutation = "_".join([mutant_fow, mutant_rev])
                        if not mutation in dbarcode[full_barcode]:
                            dbarcode[full_barcode][mutation] = 0
                        dbarcode[full_barcode][mutation] += 1
                        count_cases['Expect_MONOmut_disagree'] += 1
            elif mcode_fow + mcode_rev in ["ON","NO"]:
                # Case 3: One side was seen WT while the other was seen mutated as expected
                if mutant_rev == "W":
                    # 3a : Mutation is in the forward
                    if not mutant_fow in dbarcode[full_barcode]:
                        dbarcode[full_barcode][mutant_fow] = 0
                    dbarcode[full_barcode][mutant_fow] += 1
                    count_cases['WT+Expect_mut'] += 1
                elif mutant_fow == "W":
                    # 3b : Mutation is in the reverse
                    if not mutant_rev in dbarcode[full_barcode]:
                        dbarcode[full_barcode][mutant_rev] = 0
                    dbarcode[full_barcode][mutant_rev] += 1
                    count_cases['WT+Expect_mut'] += 1
                else:
                    # 3c : Mutation is in the forward (>=2 means DI or TRI). The other is bullshit
                    # We tested >=1 but std_dev of synonymous increased
                    if mcode_fow == "O" and len(mutant_fow.split('.')) >= 2:
                        if not mutant_fow in dbarcode[full_barcode]:
                            dbarcode[full_barcode][mutant_fow] = 0
                        dbarcode[full_barcode][mutant_fow] += 1
                        count_cases['Expect_DITRI+Bullshit'] += 1
                    # 3d : Mutation is in the reverse (>=2 means DI or TRI). The other is bullshit
                    # We tested >=1 but std_dev of synonymous increased
                    elif mcode_rev == "O" and len(mutant_rev.split('.')) >= 2:
                        if not mutant_rev in dbarcode[full_barcode]:
                            dbarcode[full_barcode][mutant_rev] = 0
                        dbarcode[full_barcode][mutant_rev] += 1
                        count_cases['Expect_DITRI+Bullshit'] += 1
            elif mcode_fow + mcode_rev == "NN":
                # Case 4: both reads are WT
                if mutant_fow + mutant_rev == "WW":
                    if not "WW" in dbarcode[full_barcode]:
                        dbarcode[full_barcode]["WW"] = 0
                    dbarcode[full_barcode]["WW"] += 1
                    count_cases['WT_agree'] += 1
                # Case 5: Reads do not contain expected sequence
                else:
                    count_cases['NotExpected'] += 1
                    possible_fow_mut_partial = []
                    possible_rev_mut_partial = []
                    # We try to save some mutants by detecting those DI or TRI which could be canonical
                    # but spoiled by a single perturbation.
                    if len(mutant_fow.split('_')) == 2:
                        for mut_partial in mutant_fow.split('_'):
                            if len(mut_partial.split('.')) >=2 and mut_partial in self.dmut:
                                possible_fow_mut_partial.append(mut_partial)
                    if len(mutant_rev.split('_')) == 2:
                        for mut_partial in mutant_rev.split('_'):
                            if len(mut_partial.split('.')) >=2 and mut_partial in self.dmut:
                                possible_rev_mut_partial.append(mut_partial)
                    mutant_select = None
                    if len(possible_rev_mut_partial) == 1 and len(possible_fow_mut_partial)==1:
                        if possible_rev_mut_partial[0] == possible_fow_mut_partial[0]:
                            mutant_select = possible_rev_mut_partial[0]
                    if len(possible_rev_mut_partial) == 1 and len(possible_fow_mut_partial)==0:
                        mutant_select = possible_rev_mut_partial[0]
                    if len(possible_rev_mut_partial) == 0 and len(possible_fow_mut_partial)==1:
                        mutant_select = possible_fow_mut_partial[0]
                    if mutant_select:
                        if not mutant_select in dbarcode[full_barcode]:
                            dbarcode[full_barcode][mutant_select] = 0
                        dbarcode[full_barcode][mutant_select] += 1
                        count_cases['NotExpect_ContainDIorTRI_retrieved'] += 1
                    else:
                        count_cases['NotExpect_ContainDIorTRI_notkept'] += 1

                    # We keep the single mut MONO that are not in self.dmut for statistical analysis
                    if len(mutant_fow.split('_')) == 1 and len(mutant_rev.split('_')) == 1:
                        if (len(mutant_fow.split('.')) in [1,2,3] or mutant_fow == "W") and (len(mutant_rev.split('.')) in [1,2,3] or mutant_rev == "W"):
                            if mutant_fow == mutant_rev or mutant_fow == "W" or mutant_rev == "W":
                                if mutant_fow == "W":
                                    mutant_select = mutant_rev
                                elif mutant_rev == "W":
                                    mutant_select = mutant_fow
                                else:
                                    mutant_select = mutant_fow
                                if not mutant_select in dbarcode[full_barcode]:
                                    dbarcode[full_barcode][mutant_select] = 0
                                dbarcode[full_barcode][mutant_select] += 1
                                count_cases['NotExpect_{}_kept'.format(dtranslate[len(mutant_select.split('.'))])] += 1
            if len(dbarcode[full_barcode]) == 0:
                del dbarcode[full_barcode]

        ### BLOCK MPI
        if self.use_mpi:
            self.comm.Barrier()
            dbarcode = self.comm.gather(dbarcode, root=0)
            count_cases = self.comm.gather(count_cases, root=0)
            number_of_common_hexindexes = self.comm.gather(number_of_common_hexindexes, root=0)
            number_of_unmatched_hexindexes = self.comm.gather(number_of_unmatched_hexindexes, root=0)

            if self.rank == 0:
                number_of_common_hexindexes = sum(number_of_common_hexindexes)
                number_of_unmatched_hexindexes = sum(number_of_unmatched_hexindexes)
                print("> Gathered back the list of {} dictionary.".format(self.size))
                for jj, dd in enumerate(dbarcode[1:]):
                    dbarcode[0] = lt.data_merge(dbarcode[0],dd,sum_values = True)
                dbarcode = dbarcode[0]
                for dd in count_cases[1:]:
                    count_cases[0] = lt.data_merge(count_cases[0],dd,sum_values = True)
                count_cases = count_cases[0]

            dbarcode = self.comm.bcast(dbarcode, root=0)

        if self.rank == 0:
            if fileout:
                fout = open(fileout,"w")
                fout.write("# Statistics comparing forward and reverse read\n")
                fout.write("# Total number of matched pair of indexes = {}\n".format(number_of_common_hexindexes))
                fout.write("# Total number of unmatched pair of indexes = {}\n".format(number_of_unmatched_hexindexes))
                total = 0
                for elt in self.list_attributes:
                    fout.write("{} :\t{}\n".format(elt,count_cases[elt]))
                    if elt not in ['NotExpect_ContainDIorTRI_retrieved','NotExpect_ContainDIorTRI_notkept','NotExpect_MONO_kept','NotExpect_DI_kept','NotExpect_TRI_kept']:
                        total += count_cases[elt]
                fout.write("Unassigned :\t{}\n".format(number_of_common_hexindexes-total))
                fout.write("Total of assigned reads in categories :\t{}\n".format(total))

        return dbarcode

    def count_mutants(self,fileoutroot, dbarcode, code_forw_rev='assembled'):
        """
        Analyzes through the dbarcode dictionary
        Writes the barcodes in a file with their first hits.
        :return: The dictionary dmut which associates every nucmut with a list of counts
        """
        if self.rank == 0:
            print(">>> Analyzing for every mutant, the number of barcodes.")
        dmut = {}
        fout_name = '{}_{}.out'.format(fileoutroot,code_forw_rev)
        with open(fout_name,"w") as fout:
            for bc in dbarcode:
                list_mutant = []
                for mut in dbarcode[bc]:
                    list_mutant.append([dbarcode[bc][mut],mut])
                # Mutants associated to every barcodes are sorted with respect to their number of observations
                list_mutant.sort(reverse=True) # of the form [[ncount1,mut1], [ncount2,mut2], ...]
                # Barcodes for which the max number of observation is below threshold are not analyzed further
                if list_mutant[0][0] <= self.thresh_rejection_number_of_observations_per_barcode:
                    continue
                fout.write("{}\t".format(bc))
                for ii,m in enumerate(list_mutant):
                    if ii >= 4:
                        break
                    if m[1] == "W":
                        canonical = "W"
                    elif m[1] in self.dmut:
                        canonical = "O"
                    else:
                        canonical = "N"
                    fout.write("{} {} {} | ".format(m[0],canonical,m[1]))
                    # If ii == 0 We only take into account the first mutant associated with a barcode
                    # In case we have mutiple mutants per barcode
                    if ii < self.thresh_maxnumber_of_variants_per_barcode:
                        if canonical in set(["W", "O"]):
                            if ii != 0 and m[0] < self.thresh_number_of_observations_in_multimut_per_barcode:
                                # If the additional mutants are below the threshold we don't keep them
                                # Multi mutants are only considered if they are numerous
                                break
                            if not m[1] in dmut:
                                dmut[m[1]] = []
                            dmut[m[1]].append(m[0])
                fout.write("\n")
        return dmut

    def sort_barcodes_categories(self, dbarcode, fileout = None ):
        """
        Classify barcodes in different categories depending on the type of mutants and the abundance of reads associated
        Writes the barcodes in a file with their first hits.
        :return: The dictionary dbc_cat which associates a set of bacodes with every mutant categories
        """
        # Define the start and end of bc list analysis for every cpu
        size_chunks = (len(dbarcode)) / self.size + 1 #
        index_start = int(self.rank * (size_chunks))
        if self.rank == self.size - 1:
            index_stop = len(dbarcode)
        else:
            index_stop = int((self.rank + 1) * (size_chunks))


        mutcat = ["WT","MONO","DI","TRI","BAD.MONO","BAD.DI","BAD.TRI"]

        # Initialisation of the dict
        dbc_cat_ori = {key:{"{}-{}".format(ini,end):set() for ini,end in self.thresh_abundance} for key in mutcat}
        if self.rank == 0: # Very important otherwise every list won't be identical => BIG MESS
            bc_list = list(dbarcode.keys())
        else:
            bc_list = None

        ### MPI block
        if self.use_mpi:
            self.comm.Barrier()
            bc_list = self.comm.bcast(bc_list,root=0)
            list_of_dbc_cat = [dict(dbc_cat_ori) for i in range(self.size)]
            dbc_cat = self.comm.scatter(list_of_dbc_cat, root=0)
            if self.rank==0:
                print(">>> Scattered a list of {} dictionary of barcodes lists in different mutant categories".format(self.size))
        else:
            dbc_cat = dbc_cat_ori

        if self.verbose:
            print("I'm rank {} and I will work of bc_list from {} to {}".format(self.rank,index_start,index_stop))

        for bc in bc_list[index_start:index_stop]:
            list_mutant = []
            for mut in dbarcode[bc]:
                list_mutant.append([dbarcode[bc][mut],mut])
            # Mutants associated to every barcodes are sorted with respect to their number of observations
            list_mutant.sort(reverse=True)
            # Barcodes for which the max number of observation is below threshold are not analyzed further
            if list_mutant[0][0] <= self.thresh_rejection_number_of_observations_per_barcode:
                continue
            for ii,m in enumerate(list_mutant):
                #if bc == "CGCCCCTCGCXGTGAGAAGCG":
                #    print("DEBUG",ii,m,list_mutant)
                if ii >= self.thresh_maxnumber_of_variants_per_barcode:
                    break
                if m[1] == "WW":
                    for ii,threshs in enumerate(self.thresh_abundance):
                        key_thresh = "{}-{}".format(threshs[0],threshs[1])
                        if m[0] >= threshs[0] and m[0] < threshs[1]:
                            #if bc == "CGCCCCTCGCXGTGAGAAGCG":
                            #    print("DEBUG",ii,m,'WW',key_thresh)

                            dbc_cat["WT"][key_thresh].add(bc) # ii is the range of abundance
                else:
                    mut_type = '.'.join(lt.get_mutation_type(m[1],self.dmut)) # e.g. to define BAD.MONO or MONO.
                    if mut_type in mutcat:
                        for ii,threshs in enumerate(self.thresh_abundance):
                            key_thresh = "{}-{}".format(threshs[0], threshs[1])
                            if m[0] >= threshs[0] and m[0] < threshs[1]:
                                #if bc == "CGCCCCTCGCXGTGAGAAGCG":
                                #    print("DEBUG",ii,m,mut_type,key_thresh)

                                dbc_cat[mut_type][key_thresh].add(bc) # ii is the range of abundance

        ### BLOCK MPI
        if self.use_mpi:
            self.comm.Barrier()
            dbc_cat = self.comm.gather(dbc_cat, root=0)
            if self.rank == 0:
                print("> Gathered back the list of {} dictionary.".format(self.size))
                for dd in dbc_cat[1:]:
                    dbc_cat[0] = lt.data_merge(dbc_cat[0],dd)
                dbc_cat = dbc_cat[0]
            dbc_cat = self.comm.bcast(dbc_cat, root=0)

        # Print the distribution to an output file
        if self.rank == 0:
            if fileout:
                fout = open(fileout, "w")
                fout.write("# Number of barcodes for different categories of mutants found with a certain abundance\n")
                fout.write("# $1 = Categories of mutants\n")
                fout.write("# $2 = Range of read abundance first_value <= abundance < second_value\n")
                fout.write("# $3 = Numbers of barcodes\n")
                for cat in mutcat:
                    for ii,thresh in enumerate(self.thresh_abundance):
                        key_thresh = "{}-{}".format(thresh[0], thresh[1])
                        report = "{}\t{}-{}\t{}\n".format(cat,thresh[0],thresh[1], len(dbc_cat[cat][key_thresh]))
                        fout.write("{}".format(report))
                fout .close()

        return dbc_cat

    def analyze_master_mutant(self):
        """
        Generates the file in which for each mutant we retrieve the distribution of the weighted barcodes
        :return:
        """
        fout_name = self.root_fname + '_mutants_for_rev.out'
        d_obs_mutants = [self.dobsmutants]
        files = [fout_name]
        for ii in range(len(files)):
            with open(files[ii],"w") as fout:
                for mut_nuc_name in self.dmut['order']: # We thread through all the expected mutants of the experiment
                    aa_mutant_name = self.dmut[mut_nuc_name]['full_mut_AA_name'].split("@@")
                    if mut_nuc_name in d_obs_mutants[ii]:
                        distrib_largest = list(d_obs_mutants[ii][mut_nuc_name])
                        distrib_largest.sort(reverse=True)
                        str_distrib_largest = '_'.join([str(x) for x in distrib_largest[:]])
                        mutant_number = len(d_obs_mutants[ii][mut_nuc_name])
                        fout.write("{}\t{}\t{}\t{}\n".format(mut_nuc_name, aa_mutant_name[1],mutant_number, str_distrib_largest))
                    else:
                        fout.write("{}\t{}\t0\t0\n".format(mut_nuc_name, aa_mutant_name[1]))

    def print_to_file_distrib_bc_interdistance(self, dbc_cat, dir_distrib = "distrib_barcode_similarities", fileout="distrib_bcdist"):
        """
        :param dbc_cat: dictionnary of categories pointing to sets of barcodes contained in each
        categories are typically : [MONO, BAD.MONO, TRI, etc...][number_of_reads supported by the barcode]
        :param dir_distrib:
        :param fileout:
        :return:
        """
        # indexes of dbc_cat correspond to the number of reads category
        category_to_print = [] # contain 3 elements:
        # 1 - a set of bc to analyse,
        # 2 - a list of sets to compare with,
        # 3 - Boolean to assess whether the bc to analyze should be by default excluded from the refset because it is contained in the reference
        category_to_print.append(["BAD.MONO",["WT","MONO"]])
        category_to_print.append(["BAD.DI",["MONO","DI"]])
        category_to_print.append(["BAD.TRI",["DI","TRI"]])
        category_to_print.append(["MONO",["WT","MONO"]])
        category_to_print.append(["DI",["MONO","DI"]])
        category_to_print.append(["TRI",["DI","TRI"]])
        #dir_distrib = "distrib_barcode_similarities{}".format(self.tag)
        if self.use_mpi:
            self.comm.Barrier()
        if self.rank == 0 and not os.path.isdir(dir_distrib):
            os.mkdir(dir_distrib)
        for threshs in [(2, 3), (50, 10000)]:
            ref_thresh = "{}-{}".format(self.thresh_abundance[-1][0],self.thresh_abundance[-1][1])
            test_thresh = "{}-{}".format(threshs[0],threshs[1])
            if test_thresh == ref_thresh:
                exclude_bc_from_refset = True
            else:
                exclude_bc_from_refset = False
            for cat in category_to_print:
                if cat[0][0] == "B" : # if it is a BAD mutant we don't have to remove it from the ref bc provided they are not BAD
                    exclude_bc_from_refset = False
                set2test = copy.copy(dbc_cat[cat[0]][test_thresh])
                setref = copy.copy(dbc_cat[cat[1][0]][ref_thresh]) # We combine several sets if required in the ref set
                for mutcat in cat[1][1:]:
                    setref.update(dbc_cat[mutcat][ref_thresh])
                filename = os.path.join(dir_distrib, "{}_{}_{}_vs_{}_{}.out".format(fileout, cat[0],test_thresh,"x".join(cat[1]),ref_thresh))
                if self.rank == 0:
                    print(">>> Running compare_barcodes_similarities on {}_{}_vs_{}_{}".format(cat[0],test_thresh,"x".join(cat[1]),ref_thresh))
                self.compare_barcodes_similarities(set2test, setref,filename, exclude_itself = exclude_bc_from_refset)

    def compare_barcodes_similarities(self, set_barcodes2search, set_refbarcodes, fileout, exclude_itself=False):
        """
        Scans a list of barcodes and compare it to a second list of barcodes to define for each barcode of the first
        ensemble what is the minimal divergence to one of the barcode for the second

        :param set_barcodes2search: The first set of barcodes
        :param set_refbarcodes: The second set of barcodes
        :param fileout: name of the file in which solely the distribution of the minimum divergence will be dumped
        :return: Histogram of the distributions
        """

        # Define the start and end of bc list analysis for every cpu
        size_chunks = (len(set_barcodes2search)) / self.size + 1 #
        index_start = int(self.rank * (size_chunks))
        if self.rank == self.size - 1:
            index_stop = len(set_barcodes2search)
        else:
            index_stop = int((self.rank + 1) * (size_chunks))

        # Defines a unique list of ordered barcodes. Required for chunking
        if self.rank == 0:
            list_barcodes2search = list(set_barcodes2search)
        else:
            list_barcodes2search = None

        ### BLOCK MPI
        if self.use_mpi:
            self.comm.Barrier()
            list_barcodes2search = self.comm.bcast(list_barcodes2search,root=0)
            list_of_ldminid = [[] for i in range(self.size)]
            l_min_id = self.comm.scatter(list_of_ldminid, root=0)
            if self.rank==0:
                print(">>> Scattered a list of {} dictionary of {} barcodes lists in different mutant categories".format(self.size,len(l_min_id)))
        else:
            l_min_id = []

        if self.verbose:
            print("I'm rank {} and I will work on list_barcodes2search from {} to {}".format(self.rank,index_start,index_stop))

        for bc in list_barcodes2search[index_start:index_stop]:
            min_id = lt.compute_minimal_id_with_bcref(bc,set_refbarcodes, exclude_bc = exclude_itself)
            l_min_id.append(min_id)

        ### BLOCK MPI
        if self.use_mpi:
            self.comm.Barrier()
            l_min_id = self.comm.gather(l_min_id, root=0)
            if self.rank == 0:
                print(">>> Gathered back the list of {} dictionary.".format(self.size))
                for ll in l_min_id[1:]:
                    l_min_id[0] = lt.data_merge(l_min_id[0],ll)
                l_min_id = l_min_id[0]
        if self.rank == 0:
            fout = open(fileout,"w")
            fout.write("# For every inter-barcode distance category (col. 1), number of cases are reported (col. 2)\n")
            for dist in range(self.max_interbc_dist+1):
                try:
                    count_dist = l_min_id.count(dist)
                except:
                    count_dist = 0
                fout.write("{}\t{}\n".format(dist,count_dist))
            fout.close()

    def print_to_file_compare_bad_correct_mono(self, dbc_cat, dbarcode, fileout, length_protein=200):
        """
        :param dbc_cat:
        :param fileout:
        :param length_protein:
        :return:
        """
        fout2 = open(fileout ,"w")
        cat2analyze = ["MONO","BAD.MONO"]
        # in dicindex, first 0 is for aa, the three next are for nucleotide count
        dicindex = {key:[[0,0,0,0] for kk in range(length_protein)] for key in cat2analyze}
        for cat in cat2analyze:
            for bc in dbc_cat[cat]["50-10000"]:
                list_mut_with_bc = [[val,key] for key,val in dbarcode[bc].items()]
                list_mut_with_bc.sort(reverse=True)
                #print("BC_TEST",bc)
                for count,mutant in list_mut_with_bc:
                    #print("MUTANT_TEST", mutant, count)
                    len_mutaa = len(mutant.split('_'))
                    len_mutnuc = len(mutant.split('.'))
                    if len_mutnuc == 1 and len_mutaa == 1:
                        if cat == 'MONO' and mutant not in self.dmut:
                            continue
                        if cat == 'BAD.MONO' and mutant in self.dmut:
                            continue
                        list_aamut, list_nucmut = lt.get_aamut(self.dn2a, self.coding_seq.strip(), mutant)
                        if not list_aamut:
                            continue
                        indexnuc = int(mutant[:-1])%3
                        aamut = list_aamut[0]
                        indexaa = int(aamut[1:-1])
                        dicindex[cat][indexaa][0] += 1
                        dicindex[cat][indexaa][1+indexnuc] += 1
                        #if cat == 'MONO' and mutant in self.dmut:
                        #    print('MONO', mutant, list_aamut, bc, dicindex[cat][indexaa])
                        break
        for ii in range(length_protein):
            fout2.write("{}\t".format(ii))
            for cat in cat2analyze:
                index = '_'.join((str(kk) for kk in dicindex[cat][ii]))
                fout2.write("{}\t".format(index))
            fout2.write("\n")
        fout2.close()

    def compare_reverse_vs_forward(self):
        """
        Not used anymore
        :return:
        """
        dbarcode = {}
        for ii,k in enumerate(self.dreverse):
            if k in self.dforward:
                forward = self.dforward[k]
                reverse = self.dreverse[k]
                if forward[-1] == reverse[-1]: # compare the two barcodes
                    barcode = forward[-1]
                    if not barcode: # Both forward and reverse are None
                        continue
                    if not barcode in dbarcode:
                        dbarcode[barcode] = {}
                    mutation = forward[0]
                    if not mutation in dbarcode[barcode]:
                        dbarcode[barcode][mutation] = 0
                    dbarcode[barcode][mutation] += 1
            else:
                if self.verbose:
                    print("READ {} NOT FOUND in forward file".format(k))
        return dbarcode

    def dump_barcode_dictionary(self):
        """
        Not used anymore
        Writes a file listing all the barcodes.
        :return:
        """
        dmut = {}
        with open(self.fileout,"w") as fout:
            for bc in self.fdbarcode:
                list_mutant = []
                for mut in self.fdbarcode[bc]:
                    list_mutant.append([self.fdbarcode[bc][mut],mut])
                list_mutant.sort(reverse=True) # list mutant is like [[ncount1,mutname1],[ncount2,mutname2],...]
                if list_mutant[0][0] <= self.thresh_rejection_number_of_observations_per_barcode:
                    continue
                fout.write("{}\t".format(bc))
                for ii,m in enumerate(list_mutant):
                    if ii >= 4:
                        break
                    if m[1] == "W":
                        canonical = "W"
                    elif m[1] in self.dmut:
                        canonical = "O"
                    else:
                        canonical = "N"
                    fout.write("{} {} {} | ".format(m[0],canonical,m[1]))
                    if ii == 0:
                        if canonical in set(["W","O"]):
                            if not m[1] in dmut:
                                dmut[m[1]] = set()
                            dmut[m[1]].add(m[0])
                fout.write("\n")
        return dmut

    def dump_quality(self):
        thresh_low = 5
        thresh_medium = 50
        for bc in self.dbarcode:
            if bc in self.dquality:
                for mut in self.dbarcode[bc].keys():
                    if mut == "WW":
                        for qual in self.dquality[bc]["W"]:
                            if self.dbarcode[bc][mut] < thresh_low:
                                self.dquality["W"]["low"].append(qual)
                            elif self.dbarcode[bc][mut] >= thresh_low and self.dbarcode[bc][mut] < thresh_medium:
                                self.dquality["W"]["medium"].append(qual)
                            elif self.dbarcode[bc][mut] >= thresh_medium:
                                self.dquality["W"]["high"].append(qual)
                    elif mut in self.dquality[bc]:
                        mut_type = self.dqcode[len(mut.split('.'))] # S, D, T
                        for qual in self.dquality[bc][mut]:
                            if self.dbarcode[bc][mut] < thresh_low:
                                self.dquality[mut_type]["low"].append(qual)
                            elif self.dbarcode[bc][mut] >= thresh_low and self.dbarcode[bc][mut] < thresh_medium:
                                self.dquality[mut_type]["medium"].append(qual)
                            elif self.dbarcode[bc][mut] >= thresh_medium:
                                self.dquality[mut_type]["high"].append(qual)
        for category in ["W","P","S","D","T"]:
            for elt in self.dquality[category]:
                fout = open("quality_{}_{}_{}.out".format(self.root_fname, category, elt),"w")
                for qual in self.dquality[category][elt]:
                    fout.write("{}\n".format(qual))
                fout.close()


    def read_config_file(self):
        """

        :return:
        """
        config = configparser.ConfigParser(interpolation=ExtendedInterpolation())
        config.read(self.config_file)

        self.coding_seq    = lt.get_seq_from_fasta(config.get("Files", "CODING_SEQ"))
        #print("CODING_SEQUENCE:\n{}".format(self.coding_seq))
        self.amplic_seq    = lt.get_seq_from_fasta(config.get("Files", "AMPLICON_SEQ"))

        self.mutfile       = config.get("Files", "TWIST_MUT_TABLE")
        faa2nuc            = config.get("Files", "AA2NUC")
        self.da2n = pickle.load(open(faa2nuc,'rb'),encoding='bytes') # due to py2 to 3
        self.da2n = lt.convert_dic_bytes2string(self.da2n) # due to py2 to 3
        fnuc2aa   = config.get("Files", "NUC2AA")
        self.dn2a = pickle.load(open(fnuc2aa,'rb'))#,encoding='bytes')

        mcodseq2ampli = re.search(self.coding_seq, self.amplic_seq)
        self.start_codseq_in_ampliseq, self.stop_codseq_in_ampliseq = mcodseq2ampli.span()


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
    parser.add_argument('-f', '--forward_parsed', action="store", help='input forward file')
    parser.add_argument('-r', '--reverse_parsed', action="store", help='input reverse file')
    parser.add_argument('-b', '--tandem_barcode', action="store_true", help='if barcode is present in forward and reverse')
    parser.add_argument('-o', '--root_outfilename', action="store", help='output barcode analysis file')
    parser.add_argument('-q', '--quality_report', action="store_true", default=False,  help='Should a Phred quality report be generated over the different categories')
    parser.add_argument('-c', '--config_file', default=default_configfile, action="store", help='name of the configuration file')
    parser.add_argument('-l', '--log_file', default=False, action="store", help='name of the log file')
    parser.add_argument('-t', '--tag', default=None, action="store", help='string used to specify the output filenames')
    parser.add_argument('--mpi', action="store_true", default=False, help='triggers mpi running')
    parser.add_argument('--verbose', action="store_true", default=False, help="print extra info")
    parser.add_argument('--cleanup', action="store_true", default=False, help="cleanup temporary files after execution")

    if len(sys.argv) == 1:
        print("type -h or --help for more help information")
        sys.exit(1)

    args = parser.parse_args()

    if args.forward_parsed and args.reverse_parsed:
        """
        Here generate the code for generating the hhr file first
        """
        A = Assembler(args.forward_parsed,
                      args.reverse_parsed,
                      args.root_outfilename,
                      args.config_file,
                      args.quality_report,
                      args.log_file,
                      args.tag,
                      use_mpi=args.mpi,
                      tandem_barcode = args.tandem_barcode,
                      verbose = args.verbose,
                      cleanup = args.cleanup)

    else:
        print("type -h or --help for more help information")
        sys.exit(1)


if __name__ == "__main__":


    Main()
