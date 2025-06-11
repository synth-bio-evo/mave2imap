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
import pprint
from copy import deepcopy
from mpi4py import MPI
from itertools import islice

import fileinput as FI
import h5py
import warnings

warnings.filterwarnings("ignore")


try:
    import configparser
except:
    import ConfigParser as configparser
from configparser import ExtendedInterpolation



class Parser:
    def __init__(self, fastq,
                 exp_orientation,
                 region,
                 config_file,
                 use_mpi,
                 outputfile,
                 verbose = False, cleanup = True):
        self.fastq = fastq
        self.use_mpi = use_mpi
        self.exp_orientation = exp_orientation
        self.region = int(region)
        self.config_file = config_file
        self.verbose = verbose
        self.cleanup = cleanup

        self.d3to1, self.d1to3 = lt.get_dicaa()

        # Parameter definition
        self.thresh_quality_reads = 0.3
        self.use_relax_pattern_in_constant_region = False

        self.start_mpi()
        if self.rank == 0:
            print("\n>>> Counting number of lines.")
        self.number_of_lines = self.get_nlines()

        if self.rank == 0:
            print("\n>>> Reading the configuration file.")
        self.read_config_file()

        if self.rank == 0:
            print("\n>>> Mapping the index of the primers on the reference sequence.")
        self.map_indexes_inrefseq()

        if self.rank == 0:
            if self.exp_orientation == 'forward':
                print("> FORWARD PATTERN_UMI:", self.forward_umi)
                print("> FORWARD CONSTANT_PRIMER:",self.forward_cst)
            if self.exp_orientation == 'reverse':
                print("> REVERSE CONSTANT_PRIMER:",self.reverse_cst)
                print("> REVERSE PATTERN_UMI:",self.reverse_umi)

        if self.rank == 0:
            print("\n>>> Extracting the expected mutation patterns.")
        self.get_dic_expect_mutations()

        if self.rank == 0:
            print("\n>>> Parsing the fastq with {} cpus".format(self.size))
        self.split_input(outputfile, rank=self.rank, size=self.size)


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

    def get_nlines(self):
        if self.rank == 0:
            number_of_lines = lt.count_file_lines(self.fastq)
            #print("I finished")
        else:
            number_of_lines = None
        if self.use_mpi:
            self.comm.Barrier()
            number_of_lines = self.comm.bcast(number_of_lines, root=0)
        if self.rank==0:
            print("> The number of reads in the input file is: {}\n".format(number_of_lines//4))
        return number_of_lines

    def get_dic_expect_mutations(self):
        """

        :return:
        """
        self.dmut = lt.parse_twist_devis_table(self.mutfile, self.dn2a, self.d3to1, start=None, stop=None, info_codon_index=None)
        #if self.verbose:
        #   pprint.pprint(self.dmut)

    def get_position_primer(self, amplicon_seq, constant_primer, orient="forward"):
        """
        Method which reduces the primer until it matches the amplicon_seq
        :param amplicon_seq: The full amplicon sequence with its flanking regions
        :param constant_primer: The primer of interest matching parly to the amplicon sequence
        :return: The delimitations for the region of the primer which matched
        """
        primer = constant_primer
        #print(primer,amplicon_seq)
        for ii in range(len(constant_primer)):
            pat = re.compile(primer)
            match = pat.search(amplicon_seq)
            if not match:
                if orient=="forward":
                    primer = constant_primer[ii+1:]
                elif orient=="reverse":
                    primer = constant_primer[:-ii-1]
            else:
                break

        if not match:
            return None, None, None
        else:
            limits = match.span()
            if self.rank == 0:
                print ("> The constant part of the {} primer '{}' matches over the sequence '{}' in positions {}".format(
                    orient,constant_primer,primer,limits))
            return limits[0], limits[1], primer

    def map_indexes_inrefseq(self):
        """
        From the data extracted in the .ini files, this method automatically
        extracts the index where the analysis of the variants should start and stop
        :return: indexes for xxx
        """

        mcodseq2ampli = re.search(self.coding_seq, self.amplic_seq)
        try:
            self.start_codseq_in_ampliseq, self.stop_codseq_in_ampliseq = mcodseq2ampli.span()
        except AttributeError:
            print("!!! ERROR: The coding sequence and amplicon file were not properly implemented.\n"
                  "-> Amplicon should correspond to the largest DNA fragment embedding the sequences from the vector used as constant parts, etc... \n"
                  "-> Coding sequence should be entirely embedded in the amplicon sequence\n")
            sys.exit()

        start_f, stop_f, matched_f = self.get_position_primer(self.amplic_seq, self.forward_cst,orient='forward')
        start_r, stop_r, matched_r = self.get_position_primer(self.amplic_seq, self.reverse_cst,orient='reverse')
        if start_f is None:
            print("!!! ERROR: The forward primer could not be found in the amplicon. Check your inputs")
            sys.exit()
        if start_r is None:
            print("!!! ERROR: The reverse primer could not be found in the amplicon. Check your inputs")
            sys.exit()

        self.start_read_in_codseq = stop_f - self.start_codseq_in_ampliseq
        self.stop_read_in_codseq = start_r - self.start_codseq_in_ampliseq
        if self.rank==0:
            print("> Positions of Coding Sequence Matching the Amplicon Seq. : {} - {} ".format(
                self.start_codseq_in_ampliseq, self.stop_codseq_in_ampliseq))
            print("> Forward Primer in Coding Sequence Stops at : {} ".format(self.start_read_in_codseq))
            print("> Reverse Primer Starts in Coding Sequence Starts at : {} ".format(self.stop_read_in_codseq))

    def split_input(self, outputfile, rank=0, size=1):
        """Split an input with Nlines so that every rank creates its fraction and works on it
        """

        if self.exp_orientation == 'forward':
            choice_dephaser = "({})".format("|".join(self.forward_dephaser))
            pattern1 = r"{dephaser}([ATCG]{{{len_umi}}})({cst_anchor})([ATGC]+)".format(
                dephaser=choice_dephaser,
                len_umi=len(self.forward_umi) - 2,
                cst_anchor=self.forward_cst) # we remove - 2 because of barckets
            # This is a test to increase the number of reads analyzed.
            # With pattern1_relaxed, we consider only the size of the constant region and not the sequence itself.
            pattern1_relaxed = r"^{dephaser}([ATCG]{{{len_umi}}})([ATCG]{{{len_cst_anchor}}})([ATGC]+)".format(
                dephaser=choice_dephaser,
                len_umi=len(self.forward_umi) - 2,
                len_cst_anchor=len(self.forward_cst)) # we remove - 2 because of barckets
        elif self.exp_orientation == 'reverse':
            choice_dephaser = "({})".format("|".join(self.reverse_dephaser))
            pattern1 = r"([ATGC]+)({cst_anchor})([ATCG]{{{len_umi}}}){dephaser}".format(
                dephaser=choice_dephaser,
                len_umi=len(self.reverse_umi) - 2,
                cst_anchor=self.reverse_cst) # we remove - 2 because of barckets
            # This is a test to increase the number of reads analyzed.
            # With pattern1_relaxed, we consider only the size of the constant region and not the sequence itself.
            pattern1_relaxed = r"([ATGC]+)([ATCG]{{{len_cst_anchor}}})([ATCG]{{{len_umi}}}){dephaser}$".format(
                dephaser=choice_dephaser,
                len_umi=len(self.reverse_umi) - 2,
                len_cst_anchor=len(self.reverse_cst))  # we remove - 2 because of barckets
        #print(pattern1)
        if rank == 0:
            print("> Splitting the input file into {} parts".format(size))

        self.outputfile = open(outputfile + ".rank{}".format(rank), 'w')
         # Ditributing the work for every rank
        full_input = self.fastq
        size_chunks = (int(self.number_of_lines) / 4) / self.size + 1 # we want to recover blocks of 4 lines
        index_start = int(rank * (size_chunks * 4))
        if rank == size - 1:
            index_stop = None
        else:
            index_stop = int((rank + 1) * (size_chunks * 4))
        if self.verbose:
            print("!! rank {}, will parse lines from {} to {}.".format(rank, index_start, index_stop))

        hit = 0
        read_index = 0
        rank_index = 10 + rank
        if self.rank == 0:
            count = 0
        else:
            count = None
        count = self.comm.scatter([count for i in range(self.size)],root=0)
        with open(full_input) as file_obj:
            n = 1
            object_to_iterate = islice(file_obj, index_start, index_stop)
            # for line in islice(file_obj, index_start, index_stop):
            for line in object_to_iterate:
                # Parse the fastq file
                l = line[: -1]
                if l[0] == "@":
                    if self.verbose:
                        print("\n############ NEW_READ #############\n")
                    n=1
                    FLAG_REACH_END_READ = False
                    read_index += 1
                    rdint = int('%d%09d' % (rank_index, read_index))
                    rdhex = hex(rdint)
                    #print(rdhex)
                    header = l
                    n += 1
                    KEEP_READ = False
                    PRIMER_DEFECT = False
                    BARCODE_DEFECT = False
                elif n == 4:
                    quality = line.strip()
                    quality_score = lt.get_quality_score(quality)
                    if self.verbose:
                        print("quality_score = {:.3f}".format(quality_score))
                    n += 1
                elif n == 3:
                    n += 1
                elif n == 2:
                    n += 1
                    # extract umi
                    if self.exp_orientation=='forward':
                        seq = line.strip()
                    elif self.exp_orientation=='reverse':
                        seq = lt.get_reverse(lt.get_complement(l.strip()))

                if n == 5:
                    count +=1
                    match = re.search(pattern1, seq)
                    if match:  # primer is detected
                        KEEP_READ = True
                    else:
                        # Try to read the reads even in case of primer defects
                        # Test the cases in which the constant part would have been mutated
                        # Only the length of the constant stretch is analyzed
                        if self.use_relax_pattern_in_constant_region and quality_score >= self.thresh_quality_reads:
                            match = re.search(pattern1_relaxed, seq)
                            if match:
                                KEEP_READ = True
                            else:
                                KEEP_READ = False
                                PRIMER_DEFECT = True
                        else:
                            KEEP_READ = False
                            PRIMER_DEFECT = True
                    if match:
                        if self.exp_orientation == 'forward':
                            dephaser = match.group(1)
                            start_seq, stop_seq = match.span(4)
                            umi = dephaser + match.group(2)
                        if self.exp_orientation == 'reverse':
                            dephaser = match.group(4)
                            start_seq, stop_seq = match.span(1)
                            umi = match.group(3) + dephaser

                        if self.verbose:
                            print("First filter for constant primer detection passed")
                            print(start_seq, stop_seq,
                                                  '\n',seq[start_seq:start_seq+30],'...',seq[stop_seq-30:stop_seq])
                            print(self.start_read_in_codseq,self.start_read_in_codseq%3,self.start_read_in_codseq//3)
                            print("len_coding_seq:",len(self.coding_seq))
                            print("coding_seq:",self.coding_seq)
                            print("####\n")
                            print(self.coding_seq[self.start_read_in_codseq:self.start_read_in_codseq+30])
                            print(seq[:30])
                    #
                    # Get the mutation code in nucleotides
                    #
                    if KEEP_READ:
                        if self.exp_orientation == 'forward':
                            mcode = lt.get_mutation_in_read(seq[start_seq:],
                                                        self.coding_seq,
                                                        self.start_read_in_codseq,
                                                        thresh_multi_mut = self.thresh_multi_mut)
                        if self.exp_orientation == 'reverse':
                            #print("start_index = {}".format(self.stop_read_in_codseq-len(seq[:stop_seq])))
                            mcode = lt.get_mutation_in_read(seq[:stop_seq],
                                                        self.coding_seq,
                                                        self.stop_read_in_codseq-len(seq[:stop_seq]),
                                                        thresh_multi_mut = self.thresh_multi_mut)
                        if self.verbose:
                            print("## Keep:{}".format(header))
                            print("UMI:",umi)
                            print("SEQ:",seq[start_seq:start_seq+30])
                            print("mcode:",mcode)
                            # getaa will raise an error in nucleotide is N
                            print(lt.get_aamut(self.dn2a, self.coding_seq.strip(), mcode))
                            print("\n")
                    else:
                        if self.verbose:
                            print("## Discard:{}\n".format(header))
                    #

                    if KEEP_READ:
                        if mcode in self.dmut: # for instance compatible with twist or trinex mutational pattern
                            mcode = "O\t" + mcode
                        else:
                            mcode = "N\t" + mcode
                        self.outputfile.write("{}\t{:.2f}\t{}\t{}\n".format(rdhex, quality_score, mcode, umi))
                    else:
                        if PRIMER_DEFECT:
                            self.outputfile.write("{}\t{:.2f}\tD\tP\n".format(rdhex, quality_score))
                            #self.outputfile.write(seq+"\n")
                            #self.outputfile.write(quality+"\n")
                        elif BARCODE_DEFECT:
                            self.outputfile.write("{}\t{:.2f}\tD\tB\n".format(rdhex, quality_score))
                        else:
                            # Normally there should not be D D but this can be checked
                            self.outputfile.write("{}\t{:.2f}\tD\tD\n".format(rdhex, quality_score))

        if rank == 0:
            print(
                "Start of the mutant analysis at position in the coding sequence:{}".format(self.start_read_in_codseq))
            print("First codon analyzed indexed as aminoacid: {}".format(self.start_read_in_codseq // 3 + 1))
        print("RANK:{} -> {}".format(self.rank,count))
        if self.use_mpi:
            count = self.comm.gather(count,root=0)
            if self.rank == 0:
                sum_count = sum(count)
        else:
            sum_count = 0
        self.outputfile.close()
        if rank == 0:
            print("SUM_RANK: {}".format(sum_count))
            listfiles = []
            for r in range(size):
                listfiles.append(outputfile + ".rank{}".format(r))
            with open(outputfile,"w") as fout, fileinput.input(listfiles) as fin:
                for line in fin:
                    fout.write(line)
            for fname in listfiles:
                os.remove(fname)


    def read_config_file(self):
        """
        Recovers the parameters required to process the input files
        region: int (1,2,3,... or N) defining which region of the gene was sequenced
        All variable are defined
        """
        config = configparser.ConfigParser(interpolation=ExtendedInterpolation())
        config.read(self.config_file)
        if self.rank==0:
            print(" Path to the config file: {}\n Content: \n".format(self.config_file))
            print(open(self.config_file).read())
            print("                 ****  END CONFIG FILE ****                 ")

        self.forward_primer = config.get("PRIMER_FORWARD","region{}".format(self.region))
        self.reverse_primer = config.get("PRIMER_REVERSE","region{}".format(self.region))

        #
        # Parameters for barcode retrieving
        #
        # forward should such as "(HHNHHNH)ATTCCGGGTCTTCCATGGCG" ; reverse as CCGACCT(DNDDNDD)

        pat_umi_cst = re.compile(r"(\(\w*\))?([ATGC]+)(\(\w*\))?")
        try:
            f_umi_cst = pat_umi_cst.findall(self.forward_primer)
            r_cst_umi = pat_umi_cst.findall(self.reverse_primer)
            self.forward_cst = f_umi_cst[0][1]
            self.forward_umi = f_umi_cst[0][0]
            self.reverse_cst = r_cst_umi[0][1]
            self.reverse_umi = r_cst_umi[0][2]
        except IndexError:
            print("ERROR: Could not parse the primers. Please check primers format in .ini file.")
            sys.exit()
        #
        # Dephaser are short nucleotide stretches that can be added to the end of the amplicons.
        #
        ndephaser          = int(config.get("DEPHASER", "DEPHASER_MAX_LENGTH"))
        if ndephaser > 0:
            self.forward_dephaser = config.get("DEPHASER", "FORWARD_DEPHASER_SEQ").split()
            self.reverse_dephaser = config.get("DEPHASER", "REVERSE_DEPHASER_SEQ").split()
        else:
            self.forward_dephaser = None
            self.reverse_dephaser = None

        self.coding_seq    = lt.get_seq_from_fasta(config.get("Files", "CODING_SEQ"))
        self.amplic_seq    = lt.get_seq_from_fasta(config.get("Files", "AMPLICON_SEQ"))

        self.mutfile       = config.get("Files", "TWIST_MUT_TABLE")

        faa2nuc            = config.get("Files", "AA2NUC")
        self.da2n = pickle.load(open(faa2nuc,'rb'),encoding='bytes') # due to py2 to 3
        self.da2n = lt.convert_dic_bytes2string(self.da2n) # due to py2 to 3
        fnuc2aa   = config.get("Files", "NUC2AA")
        self.dn2a = pickle.load(open(fnuc2aa,'rb'))#,encoding='bytes')

        #nuc2abundance = config.get("Files", "NUC2ABUND")
        self.thresh_multi_mut = int(config.get("THRESHOLDS", "THRESH_MULTI_MUT"))








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
    parser.add_argument('-f', '--fastq', action="store", help='input fastq file')
    parser.add_argument('-o', '--outputfile', action="store", help='output mutant analysis file')
    parser.add_argument('-e', '--exp_orientation', choices = ['forward', 'reverse'], action="store", help='experiment is reverse or forward')
    parser.add_argument('-r', '--region_index', action="store", help='Which region of the protein targeted. Identified by an integer 1 to N.')
    parser.add_argument('-c', '--config_file', default=default_configfile, action="store", help='name of the configuration file')
    parser.add_argument('--mpi', action="store_true", default=False, help='triggers mpi running')
    parser.add_argument('--verbose', action="store_true", default=False, help="print extra info")
    parser.add_argument('--cleanup', action="store_true", default=False, help="cleanup temporary files after execution")

    if len(sys.argv) == 1:
        print("type -h or --help for more help information")
        sys.exit(1)

    args = parser.parse_args()

    if args.fastq:
        """
        Here generate the code for generating the hhr file first
        """
        P = Parser(args.fastq,
                   args.exp_orientation,
                   args.region_index,
                   args.config_file,
                   args.mpi,
                   args.outputfile,
                   verbose = args.verbose,
                   cleanup = args.cleanup)

    else:
        print("type -h or --help for more help information")
        sys.exit(1)


if __name__ == "__main__":


    Main()
