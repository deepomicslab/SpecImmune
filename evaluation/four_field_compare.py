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

gene_list = ['A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRB1']
# gene_list = ['C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRB1']
gene_list = ['DRB1']

def hla_to_numeric(hla_string):
    # Remove the gene prefix (all non-digit characters up to the asterisk)
    # and all colons.
    numeric_hla = re.sub(r'^[^*]*\*', '', hla_string).replace(':', '')
    return numeric_hla

def parse_truth(truth_file):
    genes=[]
    truth_dict={}
    with open(truth_file, 'r') as f:
        for idx, line in enumerate(f):
            if idx == 0:
                genes=line.strip().split(',')[1:]
                # remove empty element
                genes = list(filter(None, genes))
                genes = [gene.split("_")[0] for gene in genes]
                new_genes = []
                for gene in genes:
                    if gene not in new_genes:
                        new_genes.append(gene)
                genes = new_genes[:8]
                # print(genes)
            else:
                sample_name, res = line.strip().split(',')[0], line.strip().split(',')[1:]
                if len(res) != 2*len(genes):
                    res.append('')
                # print(sample_name, res)
                if sample_name not in truth_dict:
                    truth_dict[sample_name]={}
                for idx, gene in enumerate(genes):
                    # print(len(res), len(genes), genes)
                    truth_dict[sample_name][gene]=[res[2*idx], res[2*idx+1]]
    # print (truth_dict)
    return truth_dict

def parse_truth_from_align_all(align_dir="/mnt/d/HLAPro_backup/Nanopore_optimize/hgscv2_truth_bwa/"):
    all_truth_dict = {}
    for file in os.listdir(align_dir):
        if file.endswith("_extracted_HLA_align.txt"):
            truth_file = os.path.join(align_dir, file)
            sample_truth_dict, sample_name = parse_truth_from_align(truth_file)
            all_truth_dict[sample_name] = sample_truth_dict
    # print (all_truth_dict)
    return all_truth_dict

def parse_truth_from_align(truth_file, cov_cutoff = 0.95, ide_cutoff = 0.99):
    # truth_file = "/mnt/d/HLAPro_backup/Nanopore_optimize/pacbio_truth/upload/HG00096_extracted_HLA_align.txt"
    match_len_dict = defaultdict(dict)
    identity_dict = defaultdict(dict)
    sample_truth_dict = {}
    with open(truth_file, 'r') as f:
        for idx, line in enumerate(f):
            field = line.strip().split()
            sample_name = field[0]
            gene = field[1].split("-")[-1]
            # if gene == "C":
            #     print (field)
            hap = field[2]
            allele_name = field[3]
            match_len = int(field[4])
            identity = float(field[5])  
            if  identity < ide_cutoff:
                continue

            true_allele_len = allele_length_dict[allele_name]
            if match_len/true_allele_len < cov_cutoff:
                continue

            if hap not in match_len_dict[gene]:
                match_len_dict[gene][hap] = {}
                identity_dict[gene][hap] = {}

            if allele_name not in match_len_dict[gene][hap]:
                match_len_dict[gene][hap][allele_name] = match_len
                identity_dict[gene][hap][allele_name] = identity
            elif match_len > match_len_dict[gene][hap][allele_name]:
                match_len_dict[gene][hap][allele_name] = match_len
                identity_dict[gene][hap][allele_name] = identity
    
    ## for each gene, find the alleles with top 5% match_len and 5% identity
    for gene in match_len_dict:
        sample_truth_dict[gene] = [[], []]
        ## sort the alleles by match_len
        for hap in match_len_dict[gene]:
            top_alleles = select_top_alleles(match_len_dict[gene][hap], identity_dict[gene][hap], gene)
            hap_index = int(hap[-1]) -1
            sample_truth_dict[gene][hap_index] = top_alleles
            # print (top_5_percent, "/".join(top_5_percent))
    # print (sample_truth_dict)
    return sample_truth_dict, sample_name

def select_top_alleles(my_match_len_dict, my_identity_dict, gene, len_diff_cutoff = 0.01, ide_diff_cutoff = 0.05):
    sorted_match_len = sorted(my_match_len_dict.items(), key=lambda x: x[1], reverse=True)
    
    ## find the alleles with match len no shorter than len_diff_cutoff compared to the longest match len
    good_length_list = {}
    for allele, match_len in my_match_len_dict.items():
        if float(match_len) >= float(sorted_match_len[0][1]) * (1-len_diff_cutoff):
            good_length_list[allele] = my_identity_dict[allele]

    sorted_identity = sorted(good_length_list.items(), key=lambda x: x[1], reverse=True)

    # copy_flag = assess_gene_copy(gene_mean_len[gene], float(sorted_match_len[0][1]), float(sorted_identity[0][1]))
    # if not copy_flag:
    #     return []

    ## find the alleles with identity no shorter than ide_diff_cutoff compared to the highest identity
    good_identity_list = []
    for allele, identity in good_length_list.items():
        if float(identity) >= float(sorted_identity[0][1]) * (1-ide_diff_cutoff):
            good_identity_list.append(allele)
    
    return good_identity_list

def parse_hla_hla_input(input_file):
    genes=[]
    input_dict={}
    with open(input_file, 'r') as f:
        for idx, line in enumerate(f):
            if idx == 0:  
                continue
            if idx >= 1:  
                field = line.strip().split('\t')
                gene = field[0]

                gene = del_prefix(gene)

                if gene not in input_dict:
                    input_dict[gene] = []
                if len(field) >= 3:
                    input_dict[gene].append(field[2])
                else:
                    input_dict[gene].append('')
    return input_dict

def del_prefix(a):
    if a[:4] == "HLA-":
        a = a[4:]
    return a

def parse_simu_true(input_file):
    genes=[]
    input_dict={}
    with open(input_file, 'r') as f:
        for idx, line in enumerate(f):
            if idx == 0:  
                continue
            if idx >= 1:  
                field = line.strip().split('\t')
                gene = field[0]
                gene = del_prefix(gene)
                if gene not in input_dict:
                    input_dict[gene] = []
                if len(field) >= 3:
                    input_dict[gene].append(field[1])
                    input_dict[gene].append(field[2])
                # else:
                #     input_dict[gene].append('')
    return input_dict

def parse_spechla_input(input_file):
    genes=[]
    input_dict={}
    with open(input_file, 'r') as f:
        for idx, line in enumerate(f):
            if idx == 0 or idx == 1:  
                continue
            if idx >= 2:  
                field = line.strip().split('\t')
                gene = field[0].split("_")[1]
                if gene not in input_dict:
                    input_dict[gene] = []
                # print (input_file, field[2])
                if len(field) > 2:
                    input_dict[gene].append(field[2][:-1])
                else:
                    input_dict[gene].append('')
    return input_dict

def parse_spechla_clean_input(input_file):
    print (input_file)
    genes=[]
    gene_index = {}
    input_dict={}
    with open(input_file, 'r') as f:
        for idx, line in enumerate(f):
            if idx == 0:  
                continue
            elif idx == 1:
                field = line.strip().split('\t')
                for i in range(1, len(field)):

                    gene = field[i].split("_")[1]
                    gene_index[i] = gene
                    if gene not in input_dict:
                        input_dict[gene] = []

            if idx == 2:  

                field = line.strip().split('\t')
                # print (field)
                for i in range(1, len(field)):

                    gene = gene_index[i]
                    input_dict[gene].append(field[i])
    return input_dict

def parse_all_hla_hla_input(truth_dict):
    all_hla_la_result = {}
    for sample in truth_dict:
        # input_file = f"hla_nanopore/hla_la/{sample}.txt"  # HLA*LA
        # input_file = f"/mnt/d/HLAPro_backup/Nanopore_optimize/output0/fredhutch-hla-{sample}/hlala.like.results.txt"  # SpecHLA
        input_file = f"/mnt/d/HLAPro_backup/Nanopore_optimize/output6/fredhutch-hla-{sample}/fredhutch-hla-{sample}.HLA.type.result.txt"  # SpecLong
        print (input_file)
        ## check if input file exists use os
        if os.path.exists(input_file):
            input_dict = parse_hla_hla_input(input_file)
        else:
            print(f"File {input_file} does not exist")
            input_dict = {}
        all_hla_la_result[sample] = input_dict
    return all_hla_la_result

def parse_all_spleclong_pacbio_input(outdir="/mnt/d/HLAPro_backup/Nanopore_optimize/out_pac1/"):
    all_hla_la_result = {}
    ## for all folder in outdir
    for folder in os.listdir(outdir):
        sample = folder
        # input_file = os.path.join(outdir, folder, f"{sample}.HLA.type.result.txt")
        input_file = os.path.join(outdir, folder, f"hlala.like.results.txt")
        print (input_file)
        ## check if input file exists use os
        if os.path.exists(input_file):
            input_dict = parse_hla_hla_input(input_file)
        else:
            input_file = os.path.join(outdir, folder, f"{outdir}/{sample}/{sample}/hlala.like.results.txt")
            if os.path.exists(input_file):
                input_dict = parse_hla_hla_input(input_file)
            else:
                print(f"File {input_file} does not exist")
                input_dict = {}
        all_hla_la_result[sample] = input_dict
    return all_hla_la_result

def parse_all_spechla_input(truth_dict):
    all_spechla_result = {}
    for sample in truth_dict:
        # input_file = f"hla_nanopore/details/{sample}.hla.result.details.txt" ## 
        # input_file = f"/mnt/d/HLAPro_backup/Nanopore_optimize/output/fredhutch-hla-{sample}/hla.result.details.txt"  # SpecHLA
        input_file = f"/mnt/d/HLAPro_backup/Nanopore_optimize/output2/fredhutch-hla-{sample}/hla.new.result.txt"  # new
        input_dict = parse_spechla_clean_input(input_file)
        # input_dict = parse_spechla_input(input_file)
        all_spechla_result[sample] = input_dict
    return all_spechla_result

def parse_input(input_file):
    genes=[]
    input_dict={}
    with open(input_file, 'r') as f:
        for idx, line in enumerate(f):
            if line.strip().startswith('#'):
                continue
            if idx == 1:
                genes=line.strip().split('\t')[1:]
                # remove empty element
                genes_raw = line.strip().split('\t')[1:]
                # print(genes_raw)
                genes = [it.split('_')[1] for it in genes_raw]
                # remove duplicate
                gene_dedup = []
                for gene in genes:
                    if gene not in gene_dedup:
                        gene_dedup.append(gene)
            if line.strip().startswith('Sample'):
                continue
            sample_name, res = line.strip().split('\t')[:1], line.strip().split('\t')[1:]
            if sample_name[0] not in input_dict:
                input_dict[sample_name[0]]={}
            for idx, gene in enumerate(gene_dedup):
                # print(gene, idx, res)
                input_dict[sample_name[0]][gene]=[res[2*idx], res[2*idx+1]]
    return input_dict, genes

def compare_hla(res_hla, truth_hla, digits=4):
    # convert to numeric
    res_hla_numeric = []
    for hla in res_hla:
        if hla=='':
            res_hla_numeric.append('')
        hla_numeric=hla_to_numeric(hla)
        # print(hla, hla_numeric)
        if len(hla_numeric) > digits:
            hla_numeric=hla_numeric[:digits]
        res_hla_numeric.append(hla_numeric)
    truth_hla_numeric = []
    for hla in truth_hla:
        if hla=='':
            truth_hla_numeric.append('')
        if len(hla) > digits:
            hla=hla[:digits]
        truth_hla_numeric.append(hla)
    # print(res_hla_numeric, truth_hla_numeric)
    match_cnt=0
    # print(res_hla, res_hla_numeric, truth_hla_numeric)
    for hla in res_hla_numeric:
        if hla in truth_hla_numeric:
            match_cnt+=1
    return match_cnt

def compare(truth_dict, input_dict, res_genes):
    truth_count=0
    fail_cnt=0
    fail_cnt_dict=defaultdict(int)
    right_cnt_dict=defaultdict(int)
    for sample_name, items in input_dict.items():
        
        for gene, res in items.items():
            
            if gene in truth_dict[sample_name]:
                # print(sample_name, gene)
                
                truth_hla = truth_dict[sample_name][gene]
                res_hla = res
                if gene == "DPA1":
                    print(res_hla, truth_hla)

                # print(truth_hla, res_hla, compare_hla(res_hla, truth_hla))
                if compare_hla(res_hla, truth_hla) == 0:
                    # print("all fail:", sample_name, gene, res_hla, truth_hla)
                    fail_cnt_dict[gene]+=2
                    fail_cnt +=2
                elif compare_hla(res_hla, truth_hla) == 1:
                    # print("single fail:", sample_name, gene, res_hla, truth_hla)
                    fail_cnt_dict[gene]+=1
                    fail_cnt +=1
                    right_cnt_dict[gene]+=1
                    truth_count+=1
                elif compare_hla(res_hla, truth_hla) == 2:
                    right_cnt_dict[gene]+=2
                    truth_count+=2
                
                # truth_count += compare_hla(res_hla, truth_hla)
    sample_cnt=len(input_dict)   
    # print("sample_cnt:",sample_cnt)  
    # print(truth_dict)
    for gene, err_cnt in fail_cnt_dict.items():
        print(f"{gene}err:{err_cnt}")
        print(f"{gene}right:{right_cnt_dict[gene]}, accuracy:{right_cnt_dict[gene]/(sample_cnt*2)}")
    with open(args.output, 'w') as f:
        print(truth_count, fail_cnt)
        print(len(input_dict)*len(res_genes))
        f.write(f"Accuracy: {truth_count/(len(input_dict)*len(res_genes))}\n")
        print(f"Accuracy: {truth_count/(len(input_dict)*len(res_genes))}")
    # calculate gene accuracy
        for gene, err_cnt in right_cnt_dict.items():
            print(f"{gene}:{right_cnt_dict[gene]/(sample_cnt*2)}")

def has_intersection(list1, list2):
    flag = False
    ## list1: truth, list2: answer
    set1 = set(list1)
    set2 = set(list2)
    # intersection = set1.intersection(set2)
    # return len(intersection) > 0
    for truth in set1:
        for answer in set2:
            if truth in answer:
                return True
    return False

def convert_field(mylist, digit=8):
    for j in range(len(mylist)):
        mylist[j] = convert_field_for_allele(mylist[j], digit)
    return mylist

def align_digit_2_truth(truth, mylist):  # not using
    for j in range(len(mylist)):
        mylist[j] = convert_field_for_allele(mylist[j], truth[0].split(":")*2)
    # print (mylist, truth)
    return mylist

def compare_four(truth_dict, all_hla_la_result, gene_list, digit=8):
    gene_dict = {}
    for sample in truth_dict:
        # if sample != "FH14":
        #     continue
        # for gene in truth_dict[sample]:
        for gene in gene_list:
            if gene not in truth_dict[sample]:
                # print ("truth_dict not in ", sample, gene, truth_dict[sample].keys())
                continue
            true_list = truth_dict[sample][gene]
            if gene not in all_hla_la_result[sample]:
                # print ("all_hla_la_result not in ", sample, gene, all_hla_la_result[sample])
                continue

            if true_list[0] == [] or true_list[1] == []:
                print ("copy != 2 for ", sample, gene, true_list)
                continue
            
            hla_la_list = all_hla_la_result[sample][gene]
            if hla_la_list[1] == '' and hla_la_list[0] == '':
                print ("inferred empty", sample, gene, hla_la_list)
                continue

            if true_list[1] == '':
                true_list[1] = true_list[0]
            if true_list[0] == '':
                true_list[0] = true_list[1]
            if hla_la_list[1] == '':
                hla_la_list[1] = hla_la_list[0]
            if hla_la_list[0] == '':
                hla_la_list[0] = hla_la_list[1]

            for i in range(2):
                # print (sample, gene, true_list[i], truth_dict[sample][gene])
                ## if true_list[i] is not a list, split it
                if type(true_list[i]) is not list:
                    true_list[i] = true_list[i].split("/")
                # print (hla_la_list, hla_la_list[i])
                if re.search(";", hla_la_list[i]):
                    hla_la_list[i] = hla_la_list[i].split(";")
                elif re.search(",", hla_la_list[i]):
                    hla_la_list[i] = hla_la_list[i].split(",")
                else:
                    hla_la_list[i] = [hla_la_list[i]]
                # print ("xx", hla_la_list, hla_la_list[i])
                hla_la_list[i] = [del_prefix(x) for x in hla_la_list[i]]
                true_list[i] = [del_prefix(x) for x in true_list[i]]

                # print ("yy", true_list[i] , hla_la_list[i])
                true_list[i] = convert_field(true_list[i], digit)
                hla_la_list[i] = convert_field(hla_la_list[i], digit)

            fir = 0
            if has_intersection(true_list[0], hla_la_list[0]):
                fir += 1
            if has_intersection(true_list[1], hla_la_list[1]):
                fir += 1   

            sec = 0
            if has_intersection(true_list[0], hla_la_list[1]):
                sec += 1
            if has_intersection(true_list[1], hla_la_list[0]):
                sec += 1  

            if gene not in gene_dict:
                gene_dict[gene] = [0, 0]
            gene_dict[gene][0] += max([fir, sec])
            gene_dict[gene][1] += 2

            if max([fir, sec]) != 2:
                print (sample, gene, true_list, "<<<>>>" ,hla_la_list, max([fir, sec]))


            # print (true_list, hla_la_list, fir, sec)
        # break
    for gene, items in gene_dict.items():
        # print (gene, items)
        print (gene, items[0], items[1], round(items[0]/items[1],2))

def compare_four_old(truth_dict, all_hla_la_result_old, gene_list, digit=8):
    gene_dict = {}
    for sample in truth_dict:
        # if sample != "FH14":
        #     continue
        # for gene in truth_dict[sample]:
        for gene in gene_list:
            if gene not in truth_dict[sample]:
                # print ("truth_dict not in ", sample, gene, truth_dict[sample].keys())
                continue
            true_list = truth_dict[sample][gene]
            if gene+"_1" not in all_hla_la_result_old[sample]:
                # print ("all_hla_la_result not in ", sample, gene, all_hla_la_result[sample])
                continue
            if true_list[0] == [] or true_list[1] == []:
                print ("copy != 2 for ", sample, gene, true_list)
                continue
            
            hla_la_list = [all_hla_la_result_old[sample][gene+"_1"], all_hla_la_result_old[sample][gene+"_2"]] 
            if true_list[1] == '':
                true_list[1] = true_list[0]
            if true_list[0] == '':
                true_list[0] = true_list[1]
            if hla_la_list[1] == '':
                hla_la_list[1] = hla_la_list[0]
            if hla_la_list[0] == '':
                hla_la_list[0] = hla_la_list[1]

            for i in range(2):
                # print (sample, gene, true_list[i], truth_dict[sample][gene])
                ## if true_list[i] is not a list, split it
                if type(true_list[i]) is not list:
                    true_list[i] = true_list[i].split("/")
                # print (hla_la_list, hla_la_list[i])
                if re.search(";", hla_la_list[i]):
                    hla_la_list[i] = hla_la_list[i].split(";")
                elif re.search(",", hla_la_list[i]):
                    hla_la_list[i] = hla_la_list[i].split(",")
                else:
                    hla_la_list[i] = [hla_la_list[i]]
                # print ("xx", hla_la_list, hla_la_list[i])
                hla_la_list[i] = [del_prefix(x) for x in hla_la_list[i]]
                true_list[i] = [del_prefix(x) for x in true_list[i]]

                # print ("yy", true_list[i] , hla_la_list[i])
                true_list[i] = convert_field(true_list[i], digit)
                hla_la_list[i] = convert_field(hla_la_list[i], digit)

            fir = 0
            if has_intersection(true_list[0], hla_la_list[0]):
                fir += 1
            if has_intersection(true_list[1], hla_la_list[1]):
                fir += 1   

            sec = 0
            if has_intersection(true_list[0], hla_la_list[1]):
                sec += 1
            if has_intersection(true_list[1], hla_la_list[0]):
                sec += 1  

            if gene not in gene_dict:
                gene_dict[gene] = [0, 0]
            gene_dict[gene][0] += max([fir, sec])
            gene_dict[gene][1] += 2

            # if max([fir, sec]) != 2:
            #     print (sample, gene, true_list, "<<<>>>" ,hla_la_list, max([fir, sec]))


            # print (true_list, hla_la_list, fir, sec)
        # break
    for gene, items in gene_dict.items():
        # print (gene, items)
        print (gene, round(items[0]/items[1],2))

def count_report_allele(truth_dict, all_hla_la_result):
    count_list = []
    count_gene_array = defaultdict(list)
    for sample in truth_dict:
        # for gene in truth_dict[sample]:
        for gene in gene_list:
            # print (all_hla_la_result[sample][gene])
            count_list.append(len(all_hla_la_result[sample][gene][0]))
            count_list.append(len(all_hla_la_result[sample][gene][1]))
            count_gene_array[gene].append(len(all_hla_la_result[sample][gene][0]))
            count_gene_array[gene].append(len(all_hla_la_result[sample][gene][1]))
    print ("report allele no.", np.mean(count_list), np.median(count_list), min(count_list), max(count_list))

    for gene in count_gene_array:
        print (gene, np.mean(count_gene_array[gene]), np.median(count_gene_array[gene]), min(count_gene_array[gene]), max(count_gene_array[gene]))

def assess_sim_module(truth, infer, gene_list):
    sample_truth_dict = parse_simu_true(truth)
    sample_infer_dict = parse_hla_hla_input(infer)
    # print (sample_truth_dict)
    truth_dict, infer_dict = {}, {}
    truth_dict["test"] = sample_truth_dict
    infer_dict["test"] = sample_infer_dict
    gene_list = [del_prefix(x) for x in gene_list]
    compare_four(truth_dict, infer_dict, gene_list)


def assess_sim():
    truth = "../test/test.HLA.hap.alleles.txt"
    infer = "../test/test/test_HLA/test_HLA.HLA.type.result.txt"
    sample_truth_dict = parse_simu_true(truth)
    sample_infer_dict = parse_hla_hla_input(infer)
    # print (sample_truth_dict)
    truth_dict, infer_dict = {}, {}
    truth_dict["test"] = sample_truth_dict
    infer_dict["test"] = sample_infer_dict

    gene_list = [ 'HLA-A', 'HLA-B', 'HLA-C', 'HLA-DMA', 'HLA-DMB', 'HLA-DOA', 'HLA-DOB', 'HLA-DPA1', 'HLA-DPB1', 'HLA-DPB2', 'HLA-DQA1', 'HLA-DQB1', 'HLA-DRA', 'HLA-DRB1', 'HLA-DRB3', 'HLA-DRB4', 'HLA-DRB5', 'HLA-E', 'HLA-F', 'HLA-G', 'HLA-H', 'HLA-J', 'HLA-K', 'HLA-L', 'HLA-P', 'HLA-V', 'HLA-DQA2', 'HLA-DPA2', 'HLA-N', 'HLA-S', 'HLA-T', 'HLA-U', 'HLA-W', 'MICA', 'MICB', 'TAP1', 'TAP2', 'HFE' ]
    gene_list = [del_prefix(x) for x in gene_list]
    compare_four(truth_dict, infer_dict, gene_list)

def main():
    nano_truth = "./4_field_truth.csv"
    # truth_dict=parse_truth("hla_nanopore/4_field_truth.csv")
    # input_dict, res_genes=parse_input(args.input)
    # compare(truth_dict, input_dict, res_genes)
    # input_dict = parse_hla_hla_input("hla_nanopore/hla_la/1347-4843.txt")
    # input_dict = parse_spechla_input("hla_nanopore/details/1347-4843.hla.result.details.txt")
    # print (input_dict)

    # truth_dict=parse_truth("hla_nanopore/4_field_truth.csv")
    # all_hla_la_result = parse_all_hla_hla_input(truth_dict)
    # compare_four(truth_dict, all_hla_la_result, gene_list, 8)
    # count_report_allele(truth_dict, all_hla_la_result)

    # truth_dict=parse_truth(nano_truth)
    # all_spechla_result = parse_all_spechla_input(truth_dict)
    # compare_four(truth_dict, all_spechla_result, gene_list, 8)
    # count_report_allele(truth_dict, all_spechla_result)

    truth_dict=parse_truth(nano_truth)
    all_hla_la_result = parse_all_hla_hla_input(truth_dict)
    compare_four(truth_dict, all_hla_la_result, gene_list, 8)
    count_report_allele(truth_dict, all_hla_la_result)

def main_pacbio(gene_list, truth_dir, result_dir):
    ## remove HLA- prefix in gene_list
    gene_list = [x.split("-")[-1] for x in gene_list]
    # gene_list = ['B']

    all_truth_dict = parse_truth_from_align_all(truth_dir)
    all_hla_la_result = parse_all_spleclong_pacbio_input(result_dir)

    all_old_hlala_result = parse_hlala_pacbio()
    new_truth_dict = {}
    for sample in all_hla_la_result:
        pure_sample = sample.split(".")[0]
        new_truth_dict[sample] = all_truth_dict[pure_sample]

    # print (new_truth_dict.keys(), new_truth_dict["HG00514.1"].keys())
    # print (all_hla_la_result.keys())
    compare_four(new_truth_dict, all_hla_la_result, gene_list, 8)
    print ("------------------")
    compare_four_old(new_truth_dict, all_old_hlala_result, gene_list, 8)
    # count_report_allele(all_truth_dict, all_hla_la_result)

def parse_hlala_pacbio(file_path="HLA-LA.merge.result.txt"):
    # Initialize the dictionary to hold the parsed data
    data_dict = defaultdict(dict)
    
    with open(file_path, 'r') as tsvfile:
        reader = csv.DictReader(tsvfile, delimiter='\t')
        
        for row in reader:
            sample = row['Sample']
            # For each gene, add it to the sample's dictionary
            for gene in row:
                if gene != 'Sample':
                    data_dict[sample][gene] = row[gene]
    
    return data_dict

def cal_gene_len(db_dir):
    gene_length_dict = {}
    allele_length_dict = {}
    # for each dir in the db_dir
    for gene in os.listdir(db_dir):
        gene_dir = os.path.join(db_dir, gene)
        if os.path.isdir(gene_dir):
            for file in os.listdir(gene_dir):
                if file.endswith(".fasta"):
                    fai_file = f"{gene_dir}/{file}.fai"
                    if not os.path.exists(fai_file):
                        # print (fai_file)
                        # print ("fai file does not exist")
                        return
                    gene_length_dict[gene] = []
                    f = open(fai_file, "r")
                    for line in f:
                        length = int(line.split("\t")[1])
                        allele = line.split("\t")[0]
                        allele_length_dict[allele] = length
                        # print (allele, length)
                        gene_length_dict[gene].append(length)
                    f.close()
    gene_mean_len = {}
    for gene in gene_length_dict:
        pure_gene = del_prefix(gene)
        gene_mean_len[pure_gene] = np.mean(gene_length_dict[gene])
    # print (allele_length_dict)
    return gene_mean_len, allele_length_dict

def assess_gene_copy(mean_len, max_match, max_identity, min_mat=0.8, min_identi = 0.98):
    if max_match > mean_len * min_mat and max_identity > min_identi:
        return True
    return False


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='compare results')
    # parser.add_argument('truth', help='Input VCF file path')
    # parser.add_argument('input', help='Output VCF file path')
    # parser.add_argument('output', help='Output VCF file path')

    # args = parser.parse_args()
    # main()
    # assess_sim()

    gene_class = "HLA"
    truth_dir = "/mnt/d/HLAPro_backup/Nanopore_optimize/hgscv2_truth_bwa/"
    result_dir = "/mnt/d/HLAPro_backup/Nanopore_optimize/pacbio_hla/"

    truth_dir = "/scratch/project/cs_shuaicli/wxd/hla_pacbio_new/hifi/hgscv2_truth_bwa_zip/"
    result_dir = "/scratch/project/cs_shuaicli/wxd/hla_pacbio_new/hifi/kir_typing_out/"

    db_dir = f"../db/{gene_class}/"
    gene_list, interval_dict =  get_focus_gene(gene_class)
    gene_mean_len, allele_length_dict = cal_gene_len(db_dir)
    main_pacbio(gene_list)
    
    