from collections import defaultdict
import sys, os
import pandas as pd
import numpy as np

sys.path.insert(0, sys.path[0]+'/../../scripts/')

from bed_objects import Bed_db

def load_vdj_result_bk(raw_result, sample, cutoff=0):
    """
    raw result is like:
    sample  gene    depth   phase_set       allele_1        score_1 length_1        hap_1   allele_2        score_2 length_2        hap_2   hg38_chrom      hg38_len        variant_num     hete_variant_num
    TRAV1-1 18.3    1500604 TRAV1-1*01      100.0   274     hap1    TRAV1-1*02      100.0   269     hap2    chr14   729     2       1
    """  
    store_alleles_dict = defaultdict(dict)
    store_depth_dict = defaultdict(dict)
    phase_set_dict = defaultdict(list)
    # sample = 'test'
    with open(raw_result, "r") as f:
        # skip the header
        header = f.readline()
        for line in f:
            line = line.strip().split("\t")
            gene = line[0]
            dp = float(line[1])
            if dp < cutoff:
                continue
            allele1 = line[3]
            allele2 = line[7]
            phase_set = line[2]
            # if allele1 != allele2:
            phase_set_dict[phase_set].append([allele1, allele2])

            if sample not in store_alleles_dict:
                store_alleles_dict[sample] = {}
            if gene not in store_alleles_dict[sample]:
                store_alleles_dict[sample][gene] = []
            store_alleles_dict[sample][gene] = [[allele1], [allele2]]
    
    for phase_set, alleles in phase_set_dict.items():
        print (phase_set, len(alleles))
    return store_alleles_dict

def load_vdj_result(raw_result, sample, cutoff=0):
    """
    raw result is like:
    sample  gene    depth   phase_set       allele_1        score_1 length_1        hap_1   allele_2        score_2 length_2        hap_2   hg38_chrom      hg38_len        variant_num     hete_variant_num
    TRAV1-1 18.3    1500604 TRAV1-1*01      100.0   274     hap1    TRAV1-1*02      100.0   269     hap2    chr14   729     2       1
    """  
    store_alleles_dict = defaultdict(dict)
    store_depth_dict = defaultdict(dict)
    phase_set_dict = defaultdict(list)
    hap_phase_dict = {}
    # sample = 'test'
    with open(raw_result, "r") as f:
        # skip the header
        header = f.readline()
        for line in f:
            line = line.strip().split("\t")
            gene = line[0]
            dp = float(line[1])
            if dp < cutoff:
                continue
            allele1 = line[3]
            allele2 = line[7]
            phase_set = line[2]
            # if allele1 != allele2:
            phase_set_dict[phase_set].append([allele1, allele2])
            hap_phase_dict[(allele1, allele2)] = phase_set

            if sample not in store_alleles_dict:
                store_alleles_dict[sample] = {}
            if gene not in store_alleles_dict[sample]:
                store_alleles_dict[sample][gene] = []
            store_alleles_dict[sample][gene] = [[allele1], [allele2]]
    
    support_num = []
    for phase_set, alleles in phase_set_dict.items():
        if phase_set != 'NA':
            support_num.append(len(alleles))
            # print (phase_set, len(alleles))
    # print (hap_phase_dict)
    return store_alleles_dict, support_num,hap_phase_dict

def prepare_for_vis():
    # raw_result = '/mnt/d/HLAPro_backup/Nanopore_optimize/vdj_results2/NA19240/NA19240.IG_TR_typing_result.txt'
    # load_vdj_result(raw_result)
    bed_db = Bed_db()
    hg38_gene_info = bed_db.get_hg38_gene_interval() 
    gene_list_sorted_by_pos = list(hg38_gene_info.keys())
    focus_gene_list = gene_list_sorted_by_pos[345:362]

    for sample in ['HG00512','HG00513','HG00514']:#  ['NA19238','NA19239','NA19240']:
        raw_result = f'/mnt/d/HLAPro_backup/Nanopore_optimize/vdj_results2/{sample}/{sample}.IG_TR_typing_result.txt'
        store_alleles_dict = load_vdj_result(raw_result,sample)
        for gene in focus_gene_list:
            if gene in store_alleles_dict[sample]:
                print (sample, gene, store_alleles_dict[sample][gene])
            else:
                print (sample, gene, 'NA')

def count_block_loci_num():
    support_num_data = []
    for sample in all_trio:
        raw_result = f'/mnt/d/HLAPro_backup/Nanopore_optimize/vdj_results2/{sample}/{sample}.IG_TR_typing_result.txt'
        store_alleles_dict, support_num = load_vdj_result(raw_result,sample)
        support_num_data += support_num

    print (np.mean(support_num_data), np.median(support_num_data), min(support_num_data), max(support_num_data))
    df = pd.DataFrame(support_num_data, columns = ['No. of Loci'])
    df.to_csv("block_loci_num.csv", index=False)

def check_consistence(blocks_list):
    child_haps = list(blocks_list[-1].keys())
    child_hete_haps = [x for x in child_haps if x[0] != x[1] and blocks_list[-1][x] != 'NA']
    # print (child_hete_haps)
    correct_all, num_all = 0, 0
    for i in range(1, len(child_hete_haps)):
        hap1 = child_hete_haps[i-1]
        hap2 = child_hete_haps[i]
        if blocks_list[-1][hap1] == blocks_list[-1][hap2]:
            # print (hap1, hap2, blocks_list[-1][hap1], blocks_list[-1][hap2])
            correct, num = if_trio_cons(hap1,hap2, blocks_list)
            correct_all += correct
            num_all += num
            # break
        # break
    print (correct_all, num_all, correct_all/num_all)
    return correct_all, num_all

def if_trio_cons(hap1,hap2, blocks_list):
    # p1 = []
    child_hap1 = (hap1[0], hap2[0])
    child_hap2 = (hap1[1], hap2[1])
    print ("child", child_hap1, child_hap2)
    p1, p2 = [], []
    p1_hap1, p1_hap2, p2_hap1, p2_hap2 = (), (), (), ()
    for haps in blocks_list[0]:
        if haps[0].split('*')[0] == hap1[0].split('*')[0]:
            # print (haps, blocks_list[0][haps])
            p1+= [haps, blocks_list[0][haps]]
        if haps[0].split('*')[0] == hap2[0].split('*')[0]:
            # print (haps, blocks_list[0][haps])
            p1+= [haps, blocks_list[0][haps]]
    if len(p1) == 4:
        if p1[1] == 'NA' or p1[3] == 'NA' or p1[1] == p1[3]:
            p1_hap1 = (p1[0][0], p1[2][0])
            p1_hap2 =  (p1[0][1], p1[2][1])
    print ('p1:', p1, p1_hap1, p1_hap2)  

    for haps in blocks_list[1]:
        if haps[0].split('*')[0] == hap1[0].split('*')[0]:
            # print (haps, blocks_list[1][haps])
            p2+= [haps, blocks_list[1][haps]]
        if haps[0].split('*')[0] == hap2[0].split('*')[0]:
            # print (haps, blocks_list[1][haps])
            p2+= [haps, blocks_list[1][haps]]
    if len(p2) == 4:
        if p2[1] == 'NA' or p2[3] == 'NA' or p2[1] == p2[3]:
            p2_hap1 = (p2[0][0], p2[2][0])
            p2_hap2 =  (p2[0][1], p2[2][1])
    print ('p2:', p2, p2_hap1, p2_hap2)   

    if len(p1_hap1) > 0 and len(p2_hap1) > 0:
        if child_hap1 in [p1_hap1, p1_hap2] and child_hap2 in [p2_hap1, p2_hap2]:
            print ("consistent", 2, 2)
            return 1, 1
        elif child_hap2 in [p1_hap1, p1_hap2] and child_hap1 in [p2_hap1, p2_hap2]:
            print ("consistent", 2, 2)
            return 2, 2
        elif child_hap1 in [p1_hap1, p1_hap2, p2_hap1, p2_hap2] or child_hap2 in [p1_hap1, p1_hap2, p2_hap1, p2_hap2]:
            print ("one consistent", 1, 2)
            return 1, 2
        else:
            print ("inconsistent")
            return 0, 2
    # elif len(p1_hap1) > 0 :
    #     if child_hap1 in [p1_hap1, p1_hap2] or child_hap2 in [p1_hap1, p1_hap2]:
    #         print ("consistent", 1, 1)
    #         return 1, 1
    #     else:
    #         return 0, 1
    # elif len(p2_hap1) > 0 :
    #     if child_hap1 in [p2_hap1, p2_hap2] or child_hap2 in [p2_hap1, p2_hap2]:
    #         print ("consistent", 1, 1)
    #         return 1, 1
    #     else:
    #         return 0, 1
    else:
        print ("cannot check")
        return 0, 0
    print ('\n')

        

if __name__ == "__main__":
    trio1 = ['HG00731','HG00732','HG00733']
    trio2 = ['HG00512','HG00513','HG00514']
    trio3 = ['NA19238','NA19239','NA19240']
    all_trio = trio1 + trio2 + trio3
    trio_list = [trio1, trio2, trio3]
    ## check how many blocks are trio-cosnsitent
    correct_all_all, num_all_all = 0, 0
    for trio in trio_list:
        blocks_list = []
        for sample in trio:
            raw_result = f'/mnt/d/HLAPro_backup/Nanopore_optimize/vdj_results2/{sample}/{sample}.IG_TR_typing_result.txt'
            store_alleles_dict, support_num,hap_phase_dict = load_vdj_result(raw_result,sample)
            blocks_list.append(hap_phase_dict)
        
        correct_all, num_all = check_consistence(blocks_list)
        correct_all_all += correct_all
        num_all_all += num_all
    print ("how many blocks trio consistent", correct_all_all, num_all_all, correct_all_all/num_all_all)
