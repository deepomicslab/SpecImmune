"""
for a file like
CD1D chr1 158178030 158186427
CD8A chr2 86784610 86808396
IGKC chr2 88857161 88857683
IGKJ5 chr2 88860568 88860605
IGKJ4 chr2 88860886 88860922
IGKJ3 chr2 88861221 88861258
IGKJ2 chr2 88861525 88861563
IGKJ1 chr2 88861886 88861923
IGKV4-1 chr2 88885397 88886153
IGKV5-2 chr2 88897232 88897784
IGKV7-3 chr2 88915081 88915378
IGKV2-4 chr2 88931666 88932380
IGKV1-5 chr2 88947301 88947957

Extract the regions that contain genes
"""

import sys
from collections import defaultdict

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
    interval_dict = defaultdict(list)
    with open(gene_file, "r") as f:
        for line in f:
            line = line.strip().split()
            chrom = line[1]
            start = int(line[2])
            end = int(line[3])
            interval_dict[chrom].append(list((start-gap, end+gap)))

    ## merge the intervals if they overlap
    for chrom in interval_dict:
        interval_dict[chrom] = merge_intervals(interval_dict[chrom])
    return (interval_dict)

def output_bed(interval_dict, out_file):
    with open(out_file, "w") as f:
        for chrom in interval_dict:
            for interval in interval_dict[chrom]:
                f.write(f"{chrom}:{interval[0]}-{interval[1]}\n")

if __name__ == "__main__":
    gene_file =  f"{sys.path[0]}/../gene_dist//IG_TR.gene.bed"   #sys.argv[1]
    segment_bed = f"{sys.path[0]}/../gene_dist/IG_TR.segment.bed"
    interval_dict = get_gene_interval(gene_file)
    output_bed(interval_dict, segment_bed)



