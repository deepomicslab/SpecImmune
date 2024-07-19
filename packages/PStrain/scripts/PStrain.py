#!/usr/bin/env python3

from argparse import ArgumentParser
import argparse
from my_imports import *
import merge
import pipeline_V30
import logging
from datetime import datetime
import sys
from shutil import which



def is_tool(name):
    """Check whether `name` is on PATH and marked as executable."""
    return which(name) is not None

def check_input(args, scripts_dir):
    # if options.r[-3:] == ".gz":

    if not os.path.isfile(args.cfgfile):
        print ("Error: config file is not detected.")
        sys.exit(1)

    # if is_file_zipped(args.r):
    #     print ("Error: reference file should be uncompressed.")
    #     sys.exit(1)
    if not is_tool(args.bowtie2) or not is_tool(args.bowtie2_build):
        print ("Error: bowtie2 is not installed.")
        sys.exit(1)

    if not os.path.isdir(args.bowtie2db):
        print (f"Metaphlan database {args.bowtie2db} is not detected.")
        print ("obtain the Metaphlan database...\n")
        command = f"bash {scripts_dir}/collect_metaphlan_datbase.sh -x {args.metaphlan_index} -d {args.bowtie2db} -m {args.metaphlan_version}"

        print (command)
        os.system(command)
        # elif args.metaphlan_version == 4:
        #     print (f"bash {scripts_dir}/pstrain_index_metaphlan4.sh")
        #     os.system(f"bash {scripts_dir}/pstrain_index_metaphlan4.sh")
        # else:
        #     print ("wrong value for --metaphlan_version")
        #     sys.exit(1)

    ref = args.bowtie2db + "/" + args.metaphlan_index + ".fna"
    if not os.path.isfile(ref):
        print (f"Metaphlan database reference: {ref} is not detected.")
        sys.exit(1)

    if not is_tool(args.samtools):
        print ("Error: samtools is not installed.")
        sys.exit(1)


def main():

    Usage = \
    """%(prog)s [options] -c/--conf <config file> -o/--outdir <output directory> --bowtie2db <metaphlan database> -x <metaphlan db index>

    Example: python3 PStrain.py -c config.txt -o out --bowtie2db ../mpa_vOct22_CHOCOPhlAnSGB_202403/ -x mpa_vOct22_CHOCOPhlAnSGB_202403 # Metaphlan 4
             python3 PStrain.py --metaphlan_version 3  -c config.txt -o outdir/ --bowtie2db ../mpa_v31_CHOCOPhlAn_201901/ -x mpa_v31_CHOCOPhlAn_201901 # Metaphlan 3

    The config file should follow the format:

    //
    sample: [sample1_ID]
    fq1: [forward reads fastq]
    fq2: [reverse/mate reads fastq]
    //
    sample: [sample2_ID]
    fq1: [forward reads fastq]
    fq2: [reverse/mate reads fastq]
    ...

    Help information can be found by %(prog)s -h/--help, config file format for single end reads , and additional information can be found\
    in README.MD or https://github.com/wshuai294/PStrain.
    """
    scripts_dir=sys.path[0]+'/'
    parser = ArgumentParser(description="PStrain: profile strains in shotgun metagenomic sequencing reads.",prog='python3 PStrain.py',usage=Usage, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    optional=parser._action_groups.pop()
    required=parser.add_argument_group('required arguments')
    #necessary parameter
    required.add_argument("-c", "--conf",help="The configuration file of the input samples.",dest='cfgfile',metavar='')
    required.add_argument("-o", "--outdir",help="The output directory.",dest='outdir',metavar='')
    required.add_argument("--bowtie2db",help="Path to metaphlan bowtie2db.",metavar='',\
        default=scripts_dir+'../mpa_vJun23_CHOCOPhlAnSGB_202403')
    required.add_argument("-x", "--metaphlan_index",help="metaphlan bowtie2db index.",metavar='',\
        default='mpa_vJun23_CHOCOPhlAnSGB_202403')

    #alternative parameter
    optional.add_argument("--metaphlan_version",help="Metaphlan version [3 or 4]. Used to download the corresponding database. If you build the metaphlan database yourself, just ignore this.",metavar='',\
        default=4, type=int)
    optional.add_argument("--bowtie2",help="Path to bowtie2 binary.",dest='bowtie2',metavar='',default='bowtie2')
    optional.add_argument("--bowtie2-build",help="Path to bowtie2 binary.",dest='bowtie2_build',metavar='',default='bowtie2-build')
    optional.add_argument("--samtools",help="Path to samtools binary.",dest='samtools',\
        metavar='',default='samtools')
    optional.add_argument("--metaphlan",help="Path to metaphlan.",dest='metaphlan',metavar='',\
        default='metaphlan')
    optional.add_argument("--metaphlan_output_files",help="If you have run metaphlan already,\
    give metaphlan result file in each line, the order should be the same with config file. In particular, \
    '--tax_lev s' should be added while running metaphlan.",dest='metaphlan_output_files',metavar='',\
        default='')

    optional.add_argument("-p", "--proc",help="The number of process to use for parallelizing the whole pipeline, run a sample in each process."\
        ,dest='proc',metavar='',default=1, type=int)
    optional.add_argument("-n", "--nproc",help="The number of CPUs to use for parallelizing the mapping with bowtie2."\
        ,dest='nproc',metavar='',default=1, type=int)
    optional.add_argument("-w", "--weight",help="The weight of genotype frequencies while computing loss, then the weight of\
    linked read type frequencies is 1-w. The value is between 0~1.",dest='weight',metavar='',default=0.3, type=float)
    optional.add_argument( "--lambda1",help="The weight of prior knowledge while rectifying genotype frequencies. \
    The value is between 0~1. ",dest='lambda1',metavar='',default=0.0, type=float)
    optional.add_argument( "--lambda2",help="The weight of prior estimation while rectifying second order genotype frequencies. \
    The value is between 0~1. ",dest='lambda2',metavar='',default=0.0, type=float)
    optional.add_argument("--species_dp",help="The minimum depth of species to be considered in strain profiling step (default is 5).",\
        dest='species_dp',metavar='',default=5, type=int)
    # optional.add_argument("--snp_dp",help="The minimum depth of SNVs to be considered in strain profiling step (default is 6).",\
        # dest='snp_dp',metavar='',default=6, type=int)
    optional.add_argument("--snp_ratio",help="The SNVs of which the depth are less than the specific ratio of the species's mean depth would be removed.",\
        dest='snp_dp',metavar='',default=0.45, type=float)
    optional.add_argument("--qual",help="The minimum quality score of SNVs to be considered in strain profiling step.",\
        dest='qual',metavar='',default=20, type=int)
    optional.add_argument("--similarity",help="The similarity cutoff of hierachical clustering in merge step.",\
        dest='similarity',metavar='',default=0.95, type=float)
    optional.add_argument("--elbow",help="The cutoff of elbow method while identifying strains number. \
    If the loss reduction ratio is less than the cutoff, then the strains number is determined.",\
        dest='elbow',metavar='',default=0.24, type=float)
    optional.add_argument("--gatk",help="Path to gatk binary.",dest='gatk',metavar='',\
        default=scripts_dir+'/GenomeAnalysisTK.jar')
    optional.add_argument("--picard",help="Path to picard binary.",dest='picard',metavar='',\
        default=scripts_dir+'/picard.jar')
    # optional.add_argument("--dbdir_V30",help="Path to marker gene database (default is in the ../db_V30 path).",dest='dbdir_V30',metavar='',\
    #     default=scripts_dir+'../db_V30/')
    optional.add_argument("--prior",help="The file of prior knowledge of genotype frequencies in the population.\
    Not providing this is also ok.",dest='prior',metavar='',default='None')
    parser._action_groups.append(optional)
    args = parser.parse_args()


    if len(sys.argv)==1:
        print (Usage%{'prog':sys.argv[0]})
    else:
        if not os.path.exists(args.outdir):
            os.system('mkdir '+args.outdir)
        check_input(args, scripts_dir)
        give_time = datetime.now().strftime("%Y_%m_%d_%H_%M")
        logging.basicConfig(filename = args.outdir +'/'+ 'PStrain_' + give_time + '.log',\
        format='[%(asctime)s-%(filename)s-%(levelname)s:%(message)s]', level = logging.DEBUG,filemode='w')
        logging.info("Running parameters: %s"%(' '.join(sys.argv)))
        arg_list=[args.nproc,args.species_dp,args.snp_dp,"",args.picard,args.samtools,args.bowtie2_build,args.bowtie2,'',\
            args.gatk,args.metaphlan,args.proc,args.weight,args.lambda1,args.lambda2,args.prior,args.elbow,args.qual]
        pipeline_V30.multiproc(args.outdir,args.cfgfile,arg_list,args.metaphlan_output_files, args)
        merge.species_samples(args.outdir,args.cfgfile,args.similarity)

if __name__ == "__main__":
    main()
