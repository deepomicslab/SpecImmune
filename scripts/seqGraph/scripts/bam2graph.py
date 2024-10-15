from __future__ import division
import sys
import math
import pysam
import numpy as np
from sklearn.cluster import DBSCAN
import random
import re
import os
import multiprocessing
import argparse
import time
from get_raw_bkp import getInsertSize


def calCrossReads(bam_name):
    f = open(graph, "w")
    edge_dict = {}
    bamfile = pysam.AlignmentFile(filename = bam_name, mode = 'rb')
    mean, sdev, rlen = getInsertSize(bamfile)
    rlen = int(rlen)
    insert_size = int(mean + 2*sdev) + 2 * rlen

    for read in bamfile.fetch():
        if read.is_unmapped  or read.mate_is_unmapped:
            continue
        if read.mapping_quality < min_q:
            continue
        if read.has_tag('XA'):
            continue
        if (read.reference_name == read.next_reference_name):
            continue

        ref_len = int(read.reference_name.split(':')[1].split('-')[1]) - int(read.reference_name.split(':')[1].split('-')[0])
        mate_len = int(read.next_reference_name.split(':')[1].split('-')[1]) - int(read.next_reference_name.split(':')[1].split('-')[0])

        if not (abs(read.reference_start) < insert_size or abs(ref_len - read.reference_start)< insert_size):
            continue
        if not (abs(read.next_reference_start) < insert_size or abs(mate_len - read.next_reference_start)< insert_size):
            continue

        if read.is_reverse:
            refname = read.reference_name + " -"
        else:
            refname = read.reference_name + " +"
        if not read.mate_is_reverse:
            mate_ref = read.next_reference_name + " -"
        else:
            mate_ref = read.next_reference_name + " +"


        edge = refname + " " + mate_ref

        if edge not in edge_dict:
            edge_dict[edge] = 1
        else:
            edge_dict[edge] += 1
    for chrom in chrom_copy:
        # print ("#\t", chrom, chrom_copy[chrom])
        print ("SEG\t", chrom, chrom_copy[chrom], file = f)
    for edge in edge_dict:
        if edge_dict[edge] < min_dp:
            continue
        print ("JUNC\t", edge, edge_dict[edge], file = f)
        # print (edge, edge_dict[edge])

    f.close()

def cal_copy_number():
    f = open(depth_file, 'r')
    depth_dict = {}
    all_depth = []
    for line in f:
        array = line.strip().split()
        chrom = array[0]
        pos = int(array[1])
        depth = int(array[2])
        if chrom not in depth_dict:
            depth_dict[chrom] = []
        depth_dict[chrom].append(depth)
        all_depth.append(depth)
    f.close()
    median_depth = np.median(all_depth)
    chrom_copy = {}
    for chrom in depth_dict:
        chrom_copy[chrom] = round(float(np.median(depth_dict[chrom]))/median_depth)
    return chrom_copy, median_depth


min_q = 20
min_dp_ratio = 0.3
bam_name = sys.argv[1]
graph = sys.argv[2]
depth_file = sys.argv[3]
print ("start extract edges...")
chrom_copy, median_depth = cal_copy_number()
print ("median depth:", median_depth)
min_dp = min_dp_ratio * median_depth
calCrossReads(bam_name)