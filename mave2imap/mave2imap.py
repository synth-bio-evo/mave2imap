"""
    mave2imap: This module provide means for infering and mapping interface
    hotspots based on results from MAVE (Multiplexed Assays of Variant
    Effects).

    The main function of this file corresponds to the entry point of mave2imap
    executable script
"""

# Modules import
# OS related
import os
import sys
from subprocess import run as sps_run
# Utils
from datetime import datetime
import argparse
import configparser
import warnings
warnings.filterwarnings("ignore")

# Global variables
FINAL_INIT_FN = 'final_init.ini'
fastq_files = ['f1', 'f2', 'f3', 'f4']


def init_parse(ini_f):
    """
        Recovers the parameters contained in initialization file (provide by -i or --ini argument)
    """

    init_dict = {}

    config = configparser.ConfigParser(interpolation=configparser.ExtendedInterpolation())
    config.read(ini_f)
    # print("****  START CONFIG FILE CONTENT ****\n")
    # print(open(ini_f, encoding="utf-8").read())
    # print("****  END CONFIG FILE CONTENT ****\n")

# [Parameters]
    init_dict['ncpus'] = int(config.get('Parameters', 'ncpus').split('\t')[0].split(' ')[0])
    init_dict['ompi'] = config.get('Parameters', 'OMPI_MCA_btl_vader_single_copy_mechanism').split(
        '\t')[0].split(' ')[0]
    init_dict['USE_QSUB'] = bool(config.get('Parameters', 'USE_QSUB').split('\t')[0].split(' ')[0])
    init_dict['CompFastQ'] = bool(config.get('Parameters', 'CompFastQ').split('\t')[0].split(' ')[0]
                                  )
# [FASTQ_FILES]
    init_dict['f1'] = config.get('FASTQ_FILES', 'f1').split('\t')[0].split(' ')[0]
    init_dict['f2'] = config.get('FASTQ_FILES', 'f2').split('\t')[0].split(' ')[0]
    init_dict['f3'] = config.get('FASTQ_FILES', 'f3').split('\t')[0].split(' ')[0]
    init_dict['f4'] = config.get('FASTQ_FILES', 'f4').split('\t')[0].split(' ')[0]
# [TARGET_REGION]
    init_dict['ri'] = int(config.get('TARGET_REGION', 'ri').split('\t')[0].split(' ')[0])
# [Paths]
    init_dict['TMP_DIR'] = config.get('Paths', 'TMP_DIR').split('\t')[0].split(' ')[0]
    exe_dir = os.path.dirname(__file__)  # get directory path where script is located
    config.set('Paths', 'EXE_DIR', exe_dir)
    init_dict['EXE_DIR'] = exe_dir
    init_dict['PROJECT_DIR'] = config.get('Paths', 'PROJECT_DIR').split('\t')[0].split(' ')[0]
    config.set('Paths', 'TOOLS_DIR', f"{exe_dir}/tools")
    init_dict['TOOLS_DIR'] = f"{exe_dir}/tools"
    init_dict['INI_DIR'] = config.get('Paths', 'INI_DIR').split('\t')[0].split(' ')[0]
    init_dict['QSUB_EXE'] = config.get('Paths', 'QSUB_EXE').split('\t')[0].split(' ')[0]
# [Files]
    init_dict['CODING_SEQ'] = config.get('Files', 'CODING_SEQ').split('\t')[0].split(' ')[0]
    init_dict['AMPLICON_SEQ'] = config.get('Files', 'AMPLICON_SEQ').split('\t')[0].split(' ')[0]
    init_dict['TWIST_MUT_TABLE'] = config.get('Files', 'TWIST_MUT_TABLE').split(
        '\t')[0].split(' ')[0]
    init_dict['AA2NUC'] = config.get('Files', 'AA2NUC').split('\t')[0].split(' ')[0]
    init_dict['NUC2AA'] = config.get('Files', 'NUC2AA').split('\t')[0].split(' ')[0]
    init_dict['NUC2ABUND'] = config.get('Files', 'NUC2ABUND').split('\t')[0].split(' ')[0]
# [PRIMER_FORWARD]
    init_dict['primer_for_region1'] = config.get('PRIMER_FORWARD', 'region1').split(
        '\t')[0].split(' ')[0]
    init_dict['primer_for_region2'] = config.get('PRIMER_FORWARD', 'region2').split(
        '\t')[0].split(' ')[0]
    init_dict['primer_for_region3'] = config.get('PRIMER_FORWARD', 'region2').split(
        '\t')[0].split(' ')[0]
# [PRIMER_REVERSE]
    init_dict['primer_rev_region1'] = config.get('PRIMER_REVERSE', 'region1').split(
        '\t')[0].split(' ')[0]
    init_dict['primer_rev_region2'] = config.get('PRIMER_REVERSE', 'region2').split(
        '\t')[0].split(' ')[0]
    init_dict['primer_rev_region3'] = config.get('PRIMER_REVERSE', 'region2').split(
        '\t')[0].split(' ')[0]
# [DEPHASER]
    init_dict['DEPHASER_MAX_LENGTH'] = int(config.get('DEPHASER', 'DEPHASER_MAX_LENGTH'))
    init_dict['FORWARD_DEPHASER_SEQ'] = config.get('DEPHASER', 'FORWARD_DEPHASER_SEQ').split(
        '\t')[0].split(' ')[:3]
    init_dict['REVERSE_DEPHASER_SEQ'] = config.get('DEPHASER', 'REVERSE_DEPHASER_SEQ').split(
        '\t')[0].split(' ')[:3]
# [THRESHOLDS]
    init_dict['THRESH_MULTI_MUT'] = int(config.get('THRESHOLDS', 'THRESH_MULTI_MUT').split(
        '\t')[0].split(' ')[0])
    init_dict['FILTER'] = [int(n) for n in config.get('THRESHOLDS', 'FILTER').split(
        '\t')[0].split(' ')[:3]]

    # Writing our configuration file to 'final_init.ini'
    with open(FINAL_INIT_FN, 'w', encoding='utf-8') as configfile:
        config.write(configfile)

    return init_dict


def run_command(cmd):
    '''
    Run a system command
    '''
    print(f"\n{'*' * 30}\n{datetime.now()}: Running command: '{cmd}'")
    sps_run(cmd, shell=True, check=True)


def gzip(file_path, action="extract"):
    '''
    Compress/Uncompress files using gzip/gunzip
    '''
    if action.lower() == "compress":
        print(f"Compressing file {file_path}")
        cmd = f"gzip {file_path}"
        run_command(cmd)
    else:
        print(f"Uncompressing file {file_path}")
        cmd = f"gunzip {file_path}"
        run_command(cmd)


def main():
    '''
    Main script
    '''
    # Parse arguments
    usage = "Script to parse the reads from a long read experiment.\n" + \
            "output is a dictionary encoded as a hdf5 file"
    parser = argparse.ArgumentParser(usage, prefix_chars='-+')
    parser.add_argument('-i', '--ini', action="store", help='input initialization file')
    args = parser.parse_args()

    init_dict = {}

    if args.ini:
        try:
            if os.path.isfile(args.ini):
                init_dict = init_parse(args.ini)
        except FileNotFoundError(f"File {args.ini} was not found"):
            print("Initialization file required (-i). Type -h or --help for more help information")
            sys.exit(1)
    else:
        print("Initialization file required (-i). Type -h or --help for more help information")
        sys.exit(1)

    os.environ["OMPI_MCA_btl_vader_single_copy_mechanism"] = init_dict['ompi']
    # Print the current date and time
    print(f"Run Started: {datetime.now()}\n")

    # assign variable values
    ncpus = init_dict['ncpus']
    exe_dir = init_dict['EXE_DIR']
    region_index = init_dict['ri']
    ini_file = FINAL_INIT_FN
    direction_l = ['FORWARD', 'REVERSE']
    ba_l = ['BEFORE', 'AFTER']
    thresh_filter = ",".join([str(x) for x in init_dict['FILTER']])

    # Run ReadParser and ReadAssembler with appropriate flags
    for n, f in enumerate(fastq_files):
        fastq = init_dict[f]
        if fastq.split('.')[-1] == 'gz':
            gzip(fastq, action="extract")  # extract the fastq file if required
            fastq = ".".join(fastq.split('.')[:-1])
            init_dict[f] = fastq
        even_odd = n % 2
        if even_odd == 0:
            direction = direction_l[0]
        else:
            direction = direction_l[1]
        if n in [0, 1]:
            ba = ba_l[0]
        else:
            ba = ba_l[1]

        parser_cmd = f"""time mpirun -np {ncpus} python {exe_dir}/ReadParser.py -f {fastq} -o {
            direction}_{            ba} -e {direction_l[n % 2].lower()} -r {region_index} -c {
            ini_file} --mpi"""
        run_command(parser_cmd)

        # Compress FastQ files if defined in the .ini file
        if init_dict['CompFastQ']:
            if fastq.split('.')[-1] != 'gz':
                gzip(fastq, action="compress")  # extract the fastq file if required
                fastq = f"{fastq}.gz"
                init_dict[f] = fastq

        if even_odd != 0:
            # Combine assemblies using ReadAssembler
            assembly_out_fn = f'assembly_{ba.lower()}'
            assembly_cmd = f"""time mpirun -np {ncpus} python {exe_dir}/ReadAssembler.py --mpi -f {
                direction_l[0]}_{ba} -r {direction_l[1]}_{ba} -b -o {assembly_out_fn} -c {
                ini_file} -t {ba.lower()}"""
            run_command(assembly_cmd)

    # Combine assemblies using CombineMultiExp script
    print(f"""Will combine assemblies using [THRESHOLDS] filter = {thresh_filter
                                                                 }""")
    combiner_cmd = f"""time mpirun -np {ncpus} python {exe_dir
        }/CombineMultiExp.py -b assembly_{ba_l[0].lower()} -a assembly_{ba_l[1].lower()} -o result_thresh{
        thresh_filter.replace(',', '_')} -c {ini_file} -t {thresh_filter}"""
    run_command(combiner_cmd)
    combiner_cmd = f"""time mpirun -np {ncpus} python {exe_dir}/CombineMultiExp.py -b assembly_{
        ba_l[0].lower()} -a assembly_{ba_l[1].lower()} -o result_sum_all -c {
        ini_file} -t sum_all"""
    run_command(combiner_cmd)

    # Print the current date and time
    print("\nRun ended:", datetime.now())


if __name__ == "__main__":
    main()
