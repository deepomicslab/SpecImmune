"""
Choose best-mapped allele as reference

cmd: python select_reference.py -n fredhutch-hla-FH1 -o /mnt/d/HLAPro_backup/Nanopore_optimize/output

wangshuai, wshuai294@gmail.com
"""


import os
import numpy as np
import pickle
import sys
import argparse
from collections import defaultdict


from get_allele_depth import Get_depth
from determine_gene import get_focus_gene

def run_depth(args):
    cmd = f"""
    echo search  {args["o"]}/{args["n"]}/{args["n"]}.db.bam ... 
    samtools view -bS -F 0x800  {args["o"]}/{args["n"]}/{args["n"]}.db.bam | samtools sort - >{args["o"]}/{args["n"]}/{args["n"]}.db.sort.bam
    samtools index {args["o"]}/{args["n"]}/{args["n"]}.db.sort.bam
    samtools depth -aa {args["o"]}/{args["n"]}/{args["n"]}.db.sort.bam>{args["o"]}/{args["n"]}/{args["n"]}.db.depth
    """
    os.system(cmd)
    return """%s/%s/%s.db.depth"""%(args["o"], args["n"], args["n"])

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Sort allele by depth.", add_help=False, \
    usage="python3 %(prog)s -h", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    required = parser.add_argument_group("Required arguments")
    optional = parser.add_argument_group("Optional arguments")
    # required.add_argument("-f", type=str, help="IMGT reference.", metavar="\b")
    # required.add_argument("-r", type=str, help="Long-read fastq file. PacBio or Nanopore.", metavar="\b")
    required.add_argument("-n", type=str, help="Sample ID", metavar="\b")
    required.add_argument("-o", type=str, help="The output folder to store the typing results.", metavar="\b", default="./output")
    required.add_argument("-i", type=str, help="HLA,KIR,CYP",metavar="\b", default="HLA")
    optional.add_argument("-j", type=int, help="Number of threads.", metavar="\b", default=5)
    # optional.add_argument("-m", type=int, help="Maintain this number of alleles for ILP step.", metavar="\b", default=10)
    # optional.add_argument("-b", type=float, help="The match length increase ratio higher than this value is homo [0-1].", metavar="\b", default=0.9)

    # optional.add_argument("-g", type=int, help="Whether use G group resolution annotation [0|1].", metavar="\b", default=0)
    # optional.add_argument("-m", type=int, help="1 represents typing, 0 means only read assignment", metavar="\b", default=1)
    # optional.add_argument("-y", type=str, help="Read type, [nanopore|pacbio].", metavar="\b", default="pacbio")
    # optional.add_argument("-u", type=str, help="Choose full-length or exon typing. 0 indicates full-length, 1 means exon.", metavar="\b", default="0")
    optional.add_argument("-h", "--help", action="help")
    args = vars(parser.parse_args()) 

    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(0)

    output =  """%s/%s/%s.map.txt"""%(args["o"], args["n"], args["n"])   
    gene_list, interval_dict =  get_focus_gene(args)

    depth_file = run_depth(args)
    get_depth = Get_depth(depth_file)
    get_depth.record_depth()
    record_candidate_alleles, record_allele_length = get_depth.select(output, gene_list)

    print ("result is in", output)
