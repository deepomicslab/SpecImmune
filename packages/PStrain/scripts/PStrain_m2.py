#!/usr/bin/env python3

from argparse import ArgumentParser
from my_imports import *
import merge
import pipeline
import logging
import argparse

Usage = \
"""%(prog)s [options] -c/--conf <config file> -o/--outdir <output directory>

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
required.add_argument("--dbdir",help="Path to marker gene database, can be downloaded from https://zenodo.org/doi/10.5281/zenodo.10457543.",dest='dbdir',metavar='',\
    default=scripts_dir+'../db/')
#alternative parameter
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
optional.add_argument("--species_dp",help="The minimum depth of species to be considered in strain profiling step.",\
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


optional.add_argument("--bowtie2",help="Path to bowtie2 binary.",dest='bowtie2',metavar='',default='bowtie2')
optional.add_argument("--bowtie2-build",help="Path to bowtie2 binary.",dest='bowtie2_build',metavar='',default='bowtie2-build')
optional.add_argument("--samtools",help="Path to samtools binary.",dest='samtools',\
    metavar='',default='samtools')
optional.add_argument("--metaphlan",help="Path to metaphlan.",dest='metaphlan2',metavar='',\
    default=scripts_dir+'../packages/metaphlan2/metaphlan2.py')


optional.add_argument("--gatk",help="Path to gatk binary.",dest='gatk',metavar='',\
    default=scripts_dir+'/GenomeAnalysisTK.jar')
optional.add_argument("--picard",help="Path to picard binary.",dest='picard',metavar='',\
    default=scripts_dir+'/picard.jar')
optional.add_argument("--prior",help="The file of prior knowledge of genotype frequencies in the population.",dest='prior',metavar='',\
    default=scripts_dir+'../db/prior_beta.pickle')
optional.add_argument("--skip_prior",help="Whether use prior knowledge of genotype frequencies.",metavar='',\
    default=False, type=bool)
parser._action_groups.append(optional)
args = parser.parse_args()


if len(sys.argv)==1:
    print (Usage%{'prog':sys.argv[0]})
else:
    arg_list=[args.nproc,args.species_dp,args.snp_dp,args.dbdir,args.picard,args.samtools,args.bowtie2_build,args.bowtie2,'',\
        args.gatk,args.metaphlan2,args.proc,args.weight,args.lambda1,args.lambda2,args.prior,args.elbow,args.qual]
    pipeline.multiproc(args.outdir,args.cfgfile,arg_list, args)
    merge.species_samples(args.outdir,args.cfgfile,args.similarity)
