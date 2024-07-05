import argparse
import re
from collections import defaultdict
import numpy as np
import csv

import sys, os
sys.path.insert(0, sys.path[0]+'/../scripts/')

from get_lite_db import convert_field_for_allele
from determine_gene import get_focus_gene
from db_objects import My_db
from four_field_compare import compare_four


def load_anno_result(raw_result):
    """
    raw result is like:
    sample  gene    allele  score   len     start   end     chr     hap
    IG_TR_dp50_acc98_2      IGHA1   IGHA1*02        100.000 393     1025372 1025764 chr14_igh       hap1
    IG_TR_dp50_acc98_2      IGHA2   IGHA2*02        100.000 393     905141  905533  chr14_igh       hap1
    IG_TR_dp50_acc98_2      IGHD    IGHD*01 100.000 324     1159103 1159426 chr14_igh       hap1
    IG_TR_dp50_acc98_2      IGHD2-21        IGHD2-21*02     100.000 28      1206752 1206779 chr14_igh       hap1
    IG_TR_dp50_acc98_2      IGHE    IGHE*04 100.000 330     918270  918599  chr14_igh       hap1
    """  
    store_alleles_dict = defaultdict(dict)
    sample = ''
    with open(raw_result, "r") as f:
        # skip the header
        header = f.readline()
        for line in f:
            line = line.strip().split("\t")
            sample = line[0]
            gene = line[1]
            allele = line[2]
            hap = line[-1]
            if sample not in store_alleles_dict:
                store_alleles_dict[sample] = {}
            if gene not in store_alleles_dict[sample]:
                store_alleles_dict[sample][gene] = []
            store_alleles_dict[sample][gene].append([allele])
    return store_alleles_dict


gene_list, xx =  get_focus_gene("IG_TR")

test = "/mnt/d/HLAPro_backup/Nanopore_optimize/data/sim_hap/results/IG_TR_dp50_acc98_6/IG_TR_dp50_acc98_6.IG.TR.allele.txt"
store_alleles_dict = load_anno_result(test)
# print (store_alleles_dict)
compare_four(store_alleles_dict, store_alleles_dict, gene_list)