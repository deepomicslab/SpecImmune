"""
cmd: python3 extract_VDJ_segments_from_hg38.py --hg38 /mnt/d/HLAPro_backup/Nanopore_optimize/data/hg38/GRCh38.p14.genome.fa 
"""

import sys
from collections import defaultdict
import os
import argparse

def merge_intervals(intervals):
    intervals.sort()
    merged = []
    for interval in intervals:
        if not merged or merged[-1][1] < interval[0]:
            merged.append(interval)
        else:
            merged[-1][1] = max(merged[-1][1], interval[1])
    return merged

def get_gene_interval(gene_file):
    gap = 1500000
    store_pure_gene_region = defaultdict(dict)
    interval_dict = defaultdict(list)
    with open(gene_file, "r") as f:
        for line in f:
            line = line.strip().split()
            gene = line[0]
            chrom = line[1]
            start = int(line[2])
            end = int(line[3])
            interval_dict[chrom].append(list((start-gap, end+gap)))

            store_pure_gene_region[chrom][gene] = [start, end]

    ## merge the intervals if they overlap
    for chrom in interval_dict:
        interval_dict[chrom] = merge_intervals(interval_dict[chrom])
    
    lite_gene_interval = find_gene_in_lite_ref(interval_dict, store_pure_gene_region)
    return interval_dict, lite_gene_interval

def find_gene_in_lite_ref(interval_dict, store_pure_gene_region):
    lite_gene_interval = defaultdict(list)
    for chrom in interval_dict:
        for interval in interval_dict[chrom]:
            for gene in store_pure_gene_region[chrom]:
                if store_pure_gene_region[chrom][gene][0] >= interval[0] and store_pure_gene_region[chrom][gene][1] <= interval[1]:
                    contig = f"{chrom}_{interval[0]}_{interval[1]}"
                    new_start = store_pure_gene_region[chrom][gene][0] - interval[0]
                    new_end = store_pure_gene_region[chrom][gene][1] - interval[0]
                    lite_gene_interval[gene] = [contig, new_start, new_end]
    return lite_gene_interval

def output_bed(interval_dict, out_file):
    with open(out_file, "w") as f:
        for chrom in interval_dict:
            for interval in interval_dict[chrom]:
                f.write(f"{chrom}:{interval[0]}-{interval[1]}\n")

def output_new_gene_file(lite_gene_interval, lite_gene_file):
    with open(lite_gene_file, "w") as f:
        for gene in lite_gene_interval:

            f.write(f"{gene}\t{lite_gene_interval[gene][0]}\t{lite_gene_interval[gene][1]}\t{lite_gene_interval[gene][2]}\n")

def extract_segments(hg38, segment_bed, segment):
    command = f"samtools faidx {hg38} -r {segment_bed} >{segment}"
    os.system(command)
    rename_contig(segment)
    ## index the segment with samtools and bwa
    os.system(f"samtools faidx {segment}")
    os.system(f"bwa index {segment}")

def rename_contig(segment):
    # given a segment fasta, rename its contig name, and save it to a new file
    with open(segment, "r") as f:
        lines = f.readlines()
        with open(segment, "w") as f:
            for line in lines:
                if line[0] == ">":
                    line = line.replace(":", "_")
                    line = line.replace("-", "_")
                    f.write(line)
                else:
                    f.write(line)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="get small ref for VDJ from hg38.", add_help=False, \
    usage="python3 %(prog)s -h", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    required = parser.add_argument_group("Required arguments")
    optional = parser.add_argument_group("Optional arguments")
    
    required.add_argument("--hg38", type=str, help="hg38 fasta file, used by IG_TR typing.", metavar="\b")
    optional.add_argument("--lite_ref", type=str, help="lite_ref.", metavar="\b", default=sys.path[0] + "/../VDJ_ref/IG_TR.segment.fa")
    # optional.add_argument("--lite_bed", type=str, help="save gene region in the lite ref.", metavar="\b", default=sys.path[0] + "/../db/IG_TR/IG_TR.lite.bed")

    optional.add_argument("-h", "--help", action="help")
    args = vars(parser.parse_args()) 

    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(0)

    gene_file =  f"{sys.path[0]}/../gene_dist//IG_TR.gene.bed"   #sys.argv[1]
    segment_bed = f"{sys.path[0]}/../gene_dist/IG_TR.segment.bed"
    lite_gene_file =  f"{sys.path[0]}/../gene_dist//IG_TR.gene.lite.bed" 

    # hg38 = "/mnt/d/HLAPro_backup/Nanopore_optimize/data/hg38/GRCh38.p14.genome.fa"
    # segment = "/mnt/d/HLAPro_backup/Nanopore_optimize/data/hg38/IG_TR.segment.fa"

    hg38 = args["hg38"]
    segment = args["lite_ref"]

    interval_dict, lite_gene_interval = get_gene_interval(gene_file)
    output_bed(interval_dict, segment_bed)
    output_new_gene_file(lite_gene_interval, lite_gene_file)
    extract_segments(hg38, segment_bed, segment)



