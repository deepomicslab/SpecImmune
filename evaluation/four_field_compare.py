import argparse
import re
from collections import defaultdict
import numpy as np
import csv
import pandas as pd
import json

import sys, os
sys.path.insert(0, sys.path[0]+'/../scripts/')

from get_lite_db import convert_field_for_allele
from determine_gene import get_focus_gene
from db_objects import My_db

# gene_list = ['C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRB1']
# gene_list = ['DRB1']

def hla_to_numeric(hla_string):
    # Remove the gene prefix (all non-digit characters up to the asterisk)
    # and all colons.
    numeric_hla = re.sub(r'^[^*]*\*', '', hla_string).replace(':', '')
    return numeric_hla

def get_HLA_version_conversion(HLA_allele_history = "./Allelelist_history.txt"):
    map_to_latest_version = {}
    deleted_allele_num = 0
    with open(HLA_allele_history) as f:
        for line in f:
            line = line.strip()
            if line.startswith("#") or line.startswith("HLA_ID"):
                continue
            arrs = line.split(",")
            id = arrs.pop(0)
            hla = arrs[0]
            if hla == "NA":
                print ("This allele was removed in the latest version", line)
                deleted_allele_num += 1
                # continue
            for aa in arrs:
                map_to_latest_version[aa] = hla
    print ("map_to_latest_version", len(map_to_latest_version))
    print ("deleted_allele_num", deleted_allele_num)
    # sys.exit(1)
    return (map_to_latest_version)

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

def parse_truth_from_align_all(allele_length_dict, align_dir="/mnt/d/HLAPro_backup/Nanopore_optimize/hgscv2_truth_bwa/",gene_class = "HLA", len_cutoff = 0):
    all_truth_dict = {}
    for file in os.listdir(align_dir):
        if file.endswith(f"_extracted_{gene_class}_align.txt"):
            truth_file = os.path.join(align_dir, file)
            sample_truth_dict, sample_name = parse_truth_from_align(allele_length_dict, truth_file, gene_class, len_cutoff)
            all_truth_dict[sample_name] = sample_truth_dict
    # print (all_truth_dict)
    return all_truth_dict

def parse_truth_from_align(allele_length_dict, truth_file, gene_class = "HLA", len_cutoff = 0, cov_cutoff = 0.95, ide_cutoff = 0.99):
    # truth_file = "/mnt/d/HLAPro_backup/Nanopore_optimize/pacbio_truth/upload/HG00096_extracted_HLA_align.txt"
    match_len_dict = defaultdict(dict)
    identity_dict = defaultdict(dict)
    sample_truth_dict = {}
    with open(truth_file, 'r') as f:
        for idx, line in enumerate(f):
            field = line.strip().split()
            sample_name = field[0]
            if gene_class == "HLA":
                gene = field[1].split("-")[-1]
            else:
                gene = field[1]
            # if sample_name == "HG00096" and gene == "TRAV4":
            #     print (truth_file, field)
            hap = field[2]
            allele_name = field[3]
            match_len = int(field[4])
            identity = float(field[5])  
            if identity < ide_cutoff:
                continue
            true_allele_len = allele_length_dict[allele_name]
            if match_len/true_allele_len < cov_cutoff:
                continue
            if true_allele_len < len_cutoff:
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
            top_alleles = select_top_alleles(match_len_dict[gene][hap], identity_dict[gene][hap], gene_class)
            hap_index = int(hap[-1]) -1
            sample_truth_dict[gene][hap_index] = top_alleles
            # print (top_5_percent, "/".join(top_5_percent))
    # print (sample_truth_dict)
    return sample_truth_dict, sample_name

def select_top_alleles(my_match_len_dict, my_identity_dict, gene_class, len_diff_cutoff = 0.01, ide_diff_cutoff = 0.05):
    sorted_match_len = sorted(my_match_len_dict.items(), key=lambda x: x[1], reverse=True)
    if gene_class == "CYP":
        ide_diff_cutoff = 0.00001
    
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

def parse_all_hla_hla_input(truth_dict, tool = "HLA*LA"):
    all_hla_la_result = {}
    for sample in truth_dict:
        if tool ==  "HLA*LA":
            input_file = f"/mnt/d/HLAPro_backup/Nanopore_optimize/xuedong/hla_nanopore/hla_la/{sample}.txt"  # HLA*LA
        else:
            input_file = f"/mnt/d/HLAPro_backup/Nanopore_optimize/output0/fredhutch-hla-{sample}/hlala.like.results.txt"  # SpecHLA
        # input_file = f"/mnt/d/HLAPro_backup/Nanopore_optimize/output6/fredhutch-hla-{sample}/fredhutch-hla-{sample}.HLA.type.result.txt"  # SpecLong
        print (input_file)
        ## check if input file exists use os
        if os.path.exists(input_file):
            input_dict = parse_hla_hla_input(input_file)
        else:
            print(f"File {input_file} does not exist")
            input_dict = {}
        all_hla_la_result[sample] = input_dict
    return all_hla_la_result

def parse_all_spleclong_pacbio_input(gene_class, step = 1, outdir="/mnt/d/HLAPro_backup/Nanopore_optimize/out_pac1/"):
    all_hla_la_result = {}
    sample_gene_depth_dict = {}
    if step == 1:
        suffix = ".type.result.txt"
    else:
        suffix = ".final.type.result.txt"
    
    ## for all folder in outdir
    for folder in os.listdir(outdir):
        sample = folder

        if gene_class != "IG_TR":
            # input_file = os.path.join(outdir, folder, f"{sample}.HLA.type.result.txt")
            input_file = os.path.join(outdir, folder, f"hlala.like.results.txt")
            new_input_file = os.path.join(outdir, folder, f"{sample}.{gene_class}{suffix}")
            print (input_file)
            ## check if input file exists use os
            if os.path.exists(input_file):
                input_dict = parse_hla_hla_input(input_file)
            elif os.path.exists(new_input_file):
                input_dict = parse_hla_hla_input(new_input_file)
            else:
                input_file = f"{outdir}/{sample}/{sample}/hlala.like.results.txt"
                new_input_file = os.path.join(outdir, folder, f"/{sample}/{sample}.{gene_class}{suffix}")
                if os.path.exists(input_file):
                    input_dict = parse_hla_hla_input(input_file)
                elif os.path.exists(new_input_file):
                    input_dict = parse_hla_hla_input(new_input_file)
                else:
                    print(f"File {input_file} and {new_input_file} does not exist")
                    input_dict = {}
        else:
            infer = os.path.join(outdir, f"{sample}/{sample}.IG_TR_typing_result.txt")
            # infer = os.path.join(result_dir, f"{sample}/NA18506_new/NA18506_new.IG_TR_typing_result.txt")
            if os.path.exists(infer):
                print (infer)
                input_dict, sample_infer_depth_dict = load_vdj_result(infer)
                sample_gene_depth_dict[sample] = sample_infer_depth_dict
            else:
                print(f"File {infer} does not exist")
                sys.exit(1)

        all_hla_la_result[sample] = input_dict
    if gene_class == "IG_TR":
        return all_hla_la_result, sample_gene_depth_dict
    else:
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

def convert_HLA_version(mylist, map_to_latest_version):
    for j in range(len(mylist)):
        if mylist[j] in map_to_latest_version:
            mylist[j] = map_to_latest_version[mylist[j]]
        else:
            print (f"{mylist[j]} not in map_to_latest_version")
            mylist[j]="NA"
            #sys.exit(1)
    return mylist

def align_digit_2_truth(truth, mylist):  # not using
    for j in range(len(mylist)):
        mylist[j] = convert_field_for_allele(mylist[j], truth[0].split(":")*2)
    # print (mylist, truth)
    return mylist

def fill_empty(true_list):
    if true_list[1] == '':
        true_list[1] = true_list[0]
    if true_list[0] == '':
        true_list[0] = true_list[1]
    return true_list

def check_TCR_Mutation(truth_list):
    for i in range(2):
        allele = truth_list[i][0]
        field = allele.split("_")
        if len(field) == 2:
            if field[0] == "M" or field[0] == "N":
                return True
    return False

def store_results(truth_dict, all_hla_la_result, gene_list, result_file):
    data = []
    for sample in truth_dict:

        for gene in gene_list:

            if gene not in truth_dict[sample]:
                truth_dict[sample][gene] = [[], []]

            if gene not in all_hla_la_result[sample]:
                all_hla_la_result[sample][gene] = [[], []]
            for i in range(2):
                if len(truth_dict[sample][gene][i]) == 0:
                    # print ("test", truth_dict[sample][gene][i], all_hla_la_result[sample][gene][i], "/".join(all_hla_la_result[sample][gene][i]))
                    truth_dict[sample][gene][i] = 'NA'
                    
                else:
                    truth_dict[sample][gene][i] = "/".join(truth_dict[sample][gene][i])
                if len(all_hla_la_result[sample][gene][i]) == 0:
                    all_hla_la_result[sample][gene][i] = 'NA'
                else:
                    all_hla_la_result[sample][gene][i] = "/".join(all_hla_la_result[sample][gene][i])
            data.append([sample, gene, truth_dict[sample][gene][0], truth_dict[sample][gene][1], all_hla_la_result[sample][gene][0], all_hla_la_result[sample][gene][1]])
    # print (data[1])
    ## transfer data to datafram and save in a csv file
    df = pd.DataFrame(data, columns = ['sample', 'gene', 'truth_1', 'truth_2', 'infer_1', 'infer_2'])
    df.to_csv(result_file, index=False)
    print (f"Results saved in {result_file}")
                
def compare_four(truth_dict, all_hla_la_result, gene_list, digit=8, gene_class="CYP"):
    ### input dict structure: {sample: {gene: [[a,b,c], [c,d,e]]}}
    if gene_class == "HLA":
        map_to_latest_version = get_HLA_version_conversion()
    gene_dict = {}
    for sample in truth_dict:
        # print (sample)
        if sample not in all_hla_la_result:
            print (f"{sample} not in all_hla_la_result")
            continue
        for gene in gene_list:

            if gene not in truth_dict[sample]:
                print ("gene not in truth_dict", sample, gene)
                continue
            if gene not in all_hla_la_result[sample]:
                print ("gene not in infer_dict ", sample, gene)
                continue

            true_list = truth_dict[sample][gene]
            hla_la_list = all_hla_la_result[sample][gene]

            if len(true_list) != 2:
                print ("copy != 2 for truth", sample, gene, true_list)
                continue
            if true_list[0] == [] or true_list[1] == []:
                print ("copy != 2 for ", sample, gene, true_list)
                continue
            if gene_class == "IG_TR" and check_TCR_Mutation(true_list):
                continue
            if true_list[0][0] == 'NA':
                print ("empty truth", sample, gene, true_list)
                continue

            
            
            if hla_la_list[1] == '' and hla_la_list[0] == '':
                print ("inferred empty", sample, gene, hla_la_list)
                continue

            true_list = fill_empty(true_list)
            hla_la_list = fill_empty(hla_la_list)

            for i in range(2):
                # print (sample, gene, true_list[i], truth_dict[sample][gene])
                ## if true_list[i] is not a list, split it
                if type(true_list[i]) is not list:
                    true_list[i] = true_list[i].split("/")
                # print (hla_la_list, hla_la_list[i])
                
                ### if hla_la_list[i] is a list, use join
                if type(hla_la_list[i]) is list:
                    hla_la_list[i] = ";".join(hla_la_list[i])

                if re.search(";", hla_la_list[i]):
                    hla_la_list[i] = hla_la_list[i].split(";")
                elif re.search(",", hla_la_list[i]):
                    hla_la_list[i] = hla_la_list[i].split(",")
                else:
                    hla_la_list[i] = [hla_la_list[i]]

                hla_la_list[i] = [del_prefix(x) for x in hla_la_list[i]]
                true_list[i] = [del_prefix(x) for x in true_list[i]]

                if gene_class == "HLA":  ## map to the latest IMGT version
                    true_list[i] = convert_HLA_version(true_list[i], map_to_latest_version)
                    hla_la_list[i] = convert_HLA_version(hla_la_list[i], map_to_latest_version)
                
                true_list[i] = convert_field(true_list[i], digit)
                hla_la_list[i] = convert_field(hla_la_list[i], digit)

            # fir = 0
            # if has_intersection(true_list[0], hla_la_list[0]):
            #     fir += 1
            # if has_intersection(true_list[1], hla_la_list[1]):
            #     fir += 1   

            # sec = 0
            # if has_intersection(true_list[0], hla_la_list[1]):
            #     sec += 1
            # if has_intersection(true_list[1], hla_la_list[0]):
            #     sec += 1  
            fir, sec = for_rev_compare(true_list, hla_la_list, gene_class)


            if gene not in gene_dict:
                gene_dict[gene] = [0, 0]
            gene_dict[gene][0] += max([fir, sec])
            gene_dict[gene][1] += 2

            if max([fir, sec]) != 2:
                print (sample, gene, true_list, "<<<wrong>>>" ,hla_la_list, max([fir, sec]))
            # else:
            #     print (sample, gene, true_list, "<<<correct>>>" ,hla_la_list, max([fir, sec]))
        # sys.exit(1)
    gene_accuracy_dict = cal_accuracy(gene_dict)
    
    print ("truth:")
    count_report_allele(truth_dict, gene_list)
    print ("result:")
    count_report_allele(all_hla_la_result, gene_list)
    print ("finished")

    return gene_accuracy_dict  # {gene: [gene, correct, total, accuracy]}

def for_rev_compare(true_list, hla_la_list, gene_class):
    fir = 0
    sec = 0
    if gene_class != "CYP":
        if has_intersection(true_list[0], hla_la_list[0]):
            fir += 1
        if has_intersection(true_list[1], hla_la_list[1]):
            fir += 1   

        
        if has_intersection(true_list[0], hla_la_list[1]):
            sec += 1
        if has_intersection(true_list[1], hla_la_list[0]):
            sec += 1  
    else:
        if validate_star_allele(true_list[0][0], hla_la_list[0][0]):
            fir += 1
        if validate_star_allele(true_list[1][0], hla_la_list[1][0]):
            fir += 1
        if validate_star_allele(true_list[0][0], hla_la_list[1][0]):
            sec += 1
        if validate_star_allele(true_list[1][0], hla_la_list[0][0]):
            sec += 1

    return fir, sec

def cal_accuracy(gene_dict):
    gene_accuracy_dict = {}
    total_correct, total = 0, 0
    for gene, items in gene_dict.items():
        # print (gene, items)
        total_correct += items[0]
        total += items[1]
        print (gene, items[0], items[1], round(items[0]/items[1],3))
        gene_accuracy_dict[gene] = [gene, items[0], items[1], round(items[0]/items[1],3)]
    if total == 0:
        print ("total accuracy", total_correct, total, 0)
    else:
        print ("total accuracy", total_correct, total, round(total_correct/total,3))
    return gene_accuracy_dict

def count_report_allele(all_hla_la_result, gene_list):
    count_list = []
    count_gene_array = defaultdict(list)
    total_allele_num, total_NA_num = 0, 0
    for sample in all_hla_la_result:
        # for gene in truth_dict[sample]:
        for gene in gene_list:
            # print (all_hla_la_result[sample][gene])
            for i in range(2):
                if gene not in all_hla_la_result[sample]:
                    result_allele_num = 0
                else:
                    result_list = all_hla_la_result[sample][gene][i]
                    result_allele_num  = len(result_list)
                if result_allele_num == 0:
                    total_NA_num += 1
                elif result_list[0] == "NA":
                    total_NA_num += 1
                # else:
                total_allele_num += 1

                if result_allele_num > 0:
                    count_list.append(result_allele_num)
                    count_gene_array[gene].append(result_allele_num)

    print (f"report allele no. Mean: {np.mean(count_list)}, median: {np.median(count_list)}, range:[{min(count_list)}-{max(count_list)}]")
    print (f"report allele no. NA: {total_NA_num}, total allele no.: {total_allele_num}, NA ratio: {total_NA_num/(total_NA_num+total_allele_num)}")
    # for gene in count_gene_array:
    #     print (gene, np.mean(count_gene_array[gene]), np.median(count_gene_array[gene]), min(count_gene_array[gene]), max(count_gene_array[gene]))

def load_vdj_result(raw_result, cutoff=0):
    """
    raw result is like:
    sample  gene    depth   phase_set       allele_1        score_1 length_1        hap_1   allele_2        score_2 length_2        hap_2   hg38_chrom      hg38_len        variant_num     hete_variant_num
    TRAV1-1 18.3    1500604 TRAV1-1*01      100.0   274     hap1    TRAV1-1*02      100.0   269     hap2    chr14   729     2       1
    TRAV1-2 20.0    1521603;1500604 TRAV1-2*03      100.0   266     hap1    TRAV1-2*02      100.0   176     hap2    chr14   689     4       2
    TRAV2   25.7    1590846 TRAV2*01        100.0   262     hap1    TRAV2*02        100.0   259     hap2    chr14   522     3       2
    TRAV3   25.1    1590846 TRAV3*01        100.0   285     hap1    TRAV3*01        100.0   285     hap2    chr14   608     1       0
    TRAV4   18.4    1590846 TRAV4*01        100.0   276     hap1    TRAV4*01        100.0   276     hap2    chr14   830     1       0

    """  
    store_alleles_dict = defaultdict(dict)
    store_depth_dict = defaultdict(dict)
    sample = 'test'
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
            if sample not in store_alleles_dict:
                store_alleles_dict[sample] = {}
            if gene not in store_alleles_dict[sample]:
                store_alleles_dict[sample][gene] = []
            store_alleles_dict[sample][gene] = [[allele1], [allele2]]
            store_depth_dict[sample][gene] = dp
    # print (store_alleles_dict[sample])
    return store_alleles_dict[sample], store_depth_dict[sample]

def assess_sim_module(truth, infer, gene_list, gene_class="HLA"):
    sample_truth_dict = parse_simu_true(truth)
    if gene_class != "IG_TR":
        sample_infer_dict = parse_hla_hla_input(infer)
    else:
        sample_infer_dict, sample_infer_depth_dict = load_vdj_result(infer)
    # print (sample_truth_dict)
    truth_dict, infer_dict = {}, {}
    truth_dict["test"] = sample_truth_dict
    infer_dict["test"] = sample_infer_dict
    gene_list = [del_prefix(x) for x in gene_list]
    compare_four(truth_dict, infer_dict, gene_list, 8, gene_class)

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
    compare_four(truth_dict, infer_dict, gene_list,8,gene_class="HLA")

def main():  # test in Nanopore amplicon data
    gene_list = ['A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRB1']
    nano_truth = "./4_field_truth.csv"
    # truth_dict=parse_truth("hla_nanopore/4_field_truth.csv")
    # input_dict, res_genes=parse_input(args.input)
    # compare(truth_dict, input_dict, res_genes)
    # input_dict = parse_hla_hla_input("hla_nanopore/hla_la/1347-4843.txt")
    # input_dict = parse_spechla_input("hla_nanopore/details/1347-4843.hla.result.details.txt")
    # print (input_dict)

    # truth_dict=parse_truth("hla_nanopore/4_field_truth.csv")
    # all_hla_la_result = parse_all_hla_hla_input(truth_dict)
    # compare_four(truth_dict, all_hla_la_result, gene_list, 8, "HLA")
    # count_report_allele(all_hla_la_result)

    # truth_dict=parse_truth(nano_truth)
    # all_spechla_result = parse_all_spechla_input(truth_dict)
    # compare_four(truth_dict, all_spechla_result, gene_list, 8, "HLA")
    # count_report_allele(all_spechla_result)

    truth_dict=parse_truth(nano_truth)
    # all_hla_la_result = parse_all_hla_hla_input(truth_dict, "HLA*LA")
    all_hla_la_result = parse_all_hla_hla_input(truth_dict, "SpecHLA")
    compare_four(truth_dict, all_hla_la_result, gene_list, 8, "HLA")
    # count_report_allele(all_hla_la_result)

def split_IG_TR(gene_list):
    IG_list = []
    TR_list = []
    for gene in gene_list:
        if gene.startswith("IG"):
            IG_list.append(gene)
        elif gene.startswith("TR"):
            TR_list.append(gene)
        else:
            print (gene)
            sys.exit(1)
    return IG_list, TR_list

def main_pacbio(gene_list, truth_dir, result_dir, gene_class="HLA", step = 2):
    ## remove HLA- prefix in gene_list
    if gene_class == "HLA":
        gene_list = [x.split("-")[-1] for x in gene_list]
    # gene_list = ['B']
    if gene_class != "IG_TR":
        all_truth_dict = parse_truth_from_align_all(allele_length_dict, truth_dir, gene_class)
        
    else:
        len_cutoff = 0
        all_truth_dict = parse_truth_from_align_all(allele_length_dict, truth_dir, gene_class, len_cutoff)
        IG_list, TR_list = split_IG_TR(gene_list)

    # print (all_truth_dict['NA12878'])
    # return
    all_hla_la_result = parse_all_spleclong_pacbio_input(gene_class, step, result_dir, )
    # print (all_hla_la_result)
    # return
    # all_old_hlala_result = parse_hlala_pacbio()
    new_truth_dict = {}
    for sample in all_hla_la_result:
        pure_sample = sample.split(".")[0]
        new_truth_dict[sample] = all_truth_dict[pure_sample]

    # print (new_truth_dict.keys(), new_truth_dict["HG00514.1"].keys())
    # print (all_hla_la_result.keys())
    compare_four(new_truth_dict, all_hla_la_result, gene_list, 8, gene_class)
    print ("------------------")
    # if gene_class == "IG_TR":
    #     compare_four(new_truth_dict, all_hla_la_result, IG_list, 8, gene_class)
    #     print ("------------------")
    #     compare_four(new_truth_dict, all_hla_la_result, TR_list, 8, gene_class)
    #     print ("------------------")
    
    result_file = f"{benchmark_result_dir}/{gene_class}/pacbio_{gene_class}.csv"
    store_results(new_truth_dict, all_hla_la_result, gene_list, result_file)
    # count_report_allele(all_hla_la_result)

def main_vdj_hgscv(gene_list, truth_dir, result_dir, allele_length_dict, gene_class="IG_TR",step = 2):

    # print (allele_length_dict)
    len_cutoff = 0
    all_truth_dict = parse_truth_from_align_all(allele_length_dict, truth_dir, gene_class, len_cutoff)
    IG_list, TR_list = split_IG_TR(gene_list)

    all_hla_la_result, sample_gene_depth_dict = parse_all_spleclong_pacbio_input(gene_class, step, result_dir, )

    new_truth_dict = {}
    for sample in all_hla_la_result:
        pure_sample = sample.split(".")[0]
        new_truth_dict[sample] = all_truth_dict[pure_sample]

    # print (new_truth_dict.keys(), new_truth_dict["HG00514.1"].keys())
    # print (all_hla_la_result.keys())
    spec_gene_accuracy_dict = compare_four(new_truth_dict, all_hla_la_result, gene_list, 8, gene_class)
    print ("gene number", len(spec_gene_accuracy_dict))

def main_vdj_hgscv2(gene_list, truth_dir, result_dir, allele_length_dict, sum_result_file, gene_class="IG_TR",):

    # print (allele_length_dict)
    len_cutoff = 0
    cutoff = 5
    truth_dict = parse_truth_from_align_all(allele_length_dict, truth_dir, gene_class, len_cutoff)
    IG_list, TR_list = split_IG_TR(gene_list)

    new_truth_dict = {}
    tcr_gene_list = []
    sample_list = []
    infer_dict = {}
    for sample in truth_dict:
        print (sample)
        infer = os.path.join(result_dir, f"{sample}/{sample}.IG_TR_typing_result.txt")
        if os.path.exists(infer):
            print (infer)
            sample_infer_dict, sample_infer_depth_dict = load_vdj_result(infer, cutoff)
            infer_dict[sample] = sample_infer_dict
            new_truth_dict[sample] = truth_dict[sample]
            sample_list.append(sample)
            tcr_gene_list += list(set(truth_dict[sample].keys()) & set(sample_infer_dict.keys()))

    tcr_gene_list = list(set(tcr_gene_list))
    tcr_gene_list = sorted(tcr_gene_list, key=lambda x: gene_list.index(x))
    spec_gene_accuracy_dict = compare_four(new_truth_dict, infer_dict, tcr_gene_list, 8, "IG_TR")
    print ("gene number", len(tcr_gene_list), "sample number", len(sample_list))

    shared_gene_list = []
    for gene in gene_list:
        if gene in spec_gene_accuracy_dict:
            shared_gene_list.append(gene)

    data = []
    for gene in shared_gene_list:
        data.append(spec_gene_accuracy_dict[gene] + [gene[:3], gene[:4]] )

    df = pd.DataFrame(data, columns = ['gene', 'correct', 'total', 'accuracy', 'class', 'subclass'])
    df.to_csv(sum_result_file, index=False)

    new_data = class_to_chain(data)
    df = pd.DataFrame(new_data, columns = ['chain', 'correct', 'total', 'class', 'accuracy'])
    df.to_csv(sum_result_file[:-4] + "_chain.csv", index=False)

def class_to_chain(data): ## classify vdj genes
    chain_dict = {}
    for da in data:
        chain = da[0][:4]
        main_chain = da[0][:3]
        if chain not in chain_dict:
            chain_dict[chain] = [0, 0, main_chain]
        chain_dict[chain][0] += int(da[1])
        chain_dict[chain][1] += int(da[2])
    new_data = []
    for chain in chain_dict:
        new_data.append([chain] + chain_dict[chain] + [chain_dict[chain][0]/chain_dict[chain][1]])
    return new_data



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
        # print (db_dir, gene)
        gene_dir = os.path.join(db_dir, gene)
        if os.path.isdir(gene_dir):
            for file in os.listdir(gene_dir):
                if file.endswith(".fasta"):
                    fai_file = f"{gene_dir}/{file}.fai"
                    if not os.path.exists(fai_file):
                        print (fai_file, "fai file does not exist")
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

def main_TCR(result_dir, benchmark_result_dir,gene_list,gene_class, sum_result_file, cutoff):
    truth_dict = load_TCR_truth()
    new_truth_dict = {}
    tcr_gene_list = []
    sample_list = []
    infer_dict = {}
    for sample in truth_dict:
        print (sample)
        infer = os.path.join(result_dir, f"{sample}/{sample}.IG_TR_typing_result.txt")
        # infer = os.path.join(result_dir, f"{sample}/NA18506_new/NA18506_new.IG_TR_typing_result.txt")
        if os.path.exists(infer):
            print (infer)
            sample_infer_dict, sample_infer_depth_dict = load_vdj_result(infer, cutoff)
            infer_dict[sample] = sample_infer_dict
            new_truth_dict[sample] = truth_dict[sample]
            sample_list.append(sample)

            # if sample == 'NA18506':
            #     print (sample_infer_dict['TRBJ2-5'])

            tcr_gene_list += list(set(truth_dict[sample].keys()) & set(sample_infer_dict.keys()))
    # print ("truth", new_truth_dict['NA18506']['TRBJ2-5'])
    ## sort the tcr_gene_list according to the gene list
    tcr_gene_list = list(set(tcr_gene_list))
    tcr_gene_list = sorted(tcr_gene_list, key=lambda x: gene_list.index(x))
    spec_gene_accuracy_dict = compare_four(new_truth_dict, infer_dict, tcr_gene_list, 8, "IG_TR")
    print ("gene number", len(tcr_gene_list), "sample number", len(sample_list))

    result_file = f"{benchmark_result_dir}/11samples_truth_tcr_{gene_class}.csv"
    store_results(new_truth_dict, infer_dict, tcr_gene_list, result_file)

    shared_gene_list = []
    for gene in gene_list:
        if gene in spec_gene_accuracy_dict:
            shared_gene_list.append(gene)

    data = []
    for gene in shared_gene_list:
        data.append(spec_gene_accuracy_dict[gene] + [gene[:3], gene[:4]] )

    df = pd.DataFrame(data, columns = ['gene', 'correct', 'total', 'accuracy', 'class', 'subclass'])
    df.to_csv(sum_result_file, index=False)

def load_TCR_truth():
    tcr_trut_file = "tcr_truth.csv"
    ### save the truth to a dictionary, save each sample's truth to a dictionary
    """
    ,,,gene,NA18506,NA18508,NA18507,HG02059,HG02060,HG02061,NA18956,NA18517,NA10831,HG01361,HG01175
    chr7,142301134,142301432,TRBV2,01/01,01/01,01/01,01/01,01/01,01/01,01/01,01/01,01/01,01/01,01/01
    chr7,142308753,142309048,TRBV3-1,01/01,01/01,01/01,01/01,01/01,01/01,01/01,01/01,01/01,01/01,01/01
    chr7,142313371,142313666,TRBV4-1,01/01,01/01,01/01,01/01,01/01,01/01,01/01,01/01,01/01,01/01,01/01
    """
    truth_dict = {}
    with open(tcr_trut_file, 'r') as f:
        for idx, line in enumerate(f):
            field = line.strip().split(",")
            if idx == 0:
                sample_list = field
                continue
            if field[0] == "chr7_142346196_21653":
                continue
            gene = field[3]
            for i in range(4, len(field)):
                sample = sample_list[i]
                if sample not in truth_dict:
                    truth_dict[sample] = {}
                allele = field[i].split("/")
                if len(allele) == 1:
                    # allele.append(allele[0])
                    allele = ['NA', 'NA']
                # print (allele)
                truth_dict[sample][gene] = [[allele[0]], [allele[1]]]
    return truth_dict

def get_shared_sample(truth_dict, infer_dict):
    new_truth_dict = {}
    for sample in truth_dict:
        # pure_sample = sample.split(".")[0]
        if sample in infer_dict:
            new_truth_dict[sample] = truth_dict[sample]
    return new_truth_dict

def filter_depth_sample(truth_dict, spec_depth_dict, cutoff=20):
    new_truth_dict = {}
    for sample in truth_dict:
        # pure_sample = sample.split(".")[0]
        if sample in spec_depth_dict and spec_depth_dict[sample]['CYP2D6'] >= cutoff:
            new_truth_dict[sample] = truth_dict[sample]
    return new_truth_dict

def filter_depth_sample_all(truth_dict, spec_depth_dict, cutoff=20):
    new_truth_dict = {}
    for sample in truth_dict:
        # pure_sample = sample.split(".")[0]
        if sample in spec_depth_dict:
            new_truth_dict[sample] = {}
            for gene in spec_depth_dict[sample]:
                if spec_depth_dict[sample][gene] >= cutoff:
                    new_truth_dict[sample] = truth_dict[sample][gene]
    return new_truth_dict

###### for cyp
## read a json file into a dictionary
def read_pangu_result(pangu_result):
    print (pangu_result)
    # check if the file exists
    if not os.path.exists(pangu_result):
        raise FileNotFoundError(f"{pangu_result} does not exist")
    with open(pangu_result, 'r') as f:
        pangu_result = json.load(f)
    if len(pangu_result) != 1:
        #raise FileNotFoundError(f"{pangu_result} is empty")
        return {'CYP2D6':[['NA'], ['NA']]}
    #if 'diplotype' not in pangu_result[0]:
    #    raise FileNotFoundError(f"{pangu_result} has no diplotype")
    diplotype = pangu_result[0]['diplotype']
    # print (diplotype)
    # print (pangu_result[0]['haplotypes'])
    # pangu_hap_calls = []
    # pangu_alleles = set()
    # if len(pangu_result) >= 1:
    #     for hap_dict in pangu_result[0]['haplotypes']:
    #         pangu_hap_calls.append(hap_dict['call'])
    #         for allele_dict in hap_dict['alleles']:
    #             pangu_alleles.add(allele_dict['call'])
    # # print (alleles)
    # print (pangu_alleles, pangu_hap_calls)
        
    # return pangu_result, pangu_alleles
    pure_diplotype = diplotype.split()[1]
    return get_standard_diploid(pure_diplotype)

def get_standard_diploid(pure_diplotype):
    diplotype_list = pure_diplotype.split("/")
    if len(diplotype_list) == 1:
        diplotype_list.append("NA")
    # print (pure_diplotype)
    return {'CYP2D6':[[diplotype_list[0]], [diplotype_list[1]]]}

def read_spec_result(spec_result):
    # check if the file exists
    spec_result_dict = {}
    spec_gene_depth = {}
    if not os.path.exists(spec_result):
        raise FileNotFoundError(f"{spec_result} does not exist")
    pure_diplotype = 'NA'
    with open(spec_result, 'r') as f:
        for line in f:
            field = line.strip().split("\t")
            # print (field)
            if field[1] == 'CYP2D6 Diplotype:':
                pure_diplotype = field[2]
            if field[0] != "#" and field[0] != "Locus":
                gene = field[0]
                allele = field[6]  #suballele
                read_num = int(field[4])
                spec_gene_depth[gene] = read_num
                if allele != "NA":
                    allele = allele.split(".")[0] # star allele
                    allele = "*" + allele.split('*')[1]

                if gene not in spec_result_dict:
                    spec_result_dict[gene] = []
                spec_result_dict[gene].append([allele])
    # print ("#", pure_diplotype, spec_result_dict)
    return get_standard_diploid(pure_diplotype), spec_result_dict, spec_gene_depth

def load_HPRC_CYP_truth():
    cyp_hprc_truth = 'cyp/HPRC_truth.csv'
    truth_dict = defaultdict(dict)
    ## read the csv with pd
    df = pd.read_csv(cyp_hprc_truth)
    for index, row in df.iterrows():
        sample = row['Sample']
        # ref = row['CYP2D6 Reference']
        ref = row['PacBio HiFi Call']
        field = ref.split('/')
        # print (sample, ref)
        truth_dict[sample]['CYP2D6'] = [[field[0]], [field[1]]]
    # print (truth_dict)
    return truth_dict

def load_1k_CYP_truth():
    cyp_hprc_truth = 'cyp/ont_truth_merge.csv'
    truth_dict = defaultdict(dict)
    ## read the csv with pd
    for line in open(cyp_hprc_truth):
        e_field = line.strip().split()
        sample = e_field[0]
        field = e_field[1].split('/')
        # print (sample, ref)
        truth_dict[sample]['CYP2D6'] = [[field[0]], [field[1]]]
    return truth_dict

def load_GeT_RM4():
    GeT_RM_truth = {}
    ## read the excel into dataframe
    import pandas as pd
    file = "cyp/GeT_RM_truth.csv"
    df = pd.read_csv(file)
    # print (df.columns)
    ## enumerate each row
    for index, row in df.iterrows():
        sample = row['Coriell #  https://www.coriell.org/  ']

        if sample not in GeT_RM_truth:
            GeT_RM_truth[sample] = {}
        # for enumarate the df column
        for gene in df.columns:
            if len(gene.split()) != 1:
                continue
            # print (gene)

            truth = row[gene]
            ## check if truth is nan
            if pd.isna(truth):
                continue
            if truth == 'no consensus' or truth == 'duplication':
                continue
            if len(truth.split("/")) != 2:
                continue
            # print (sample, gene, truth)
            GeT_RM_truth[sample][gene] = truth.split("/")

    return GeT_RM_truth

def validate_star_allele(a, b): ## for CYP2D6
    # a truth
    # b result
    # get star allele from suballele
    # b = b.split('.')[0]
    a = a.replace(" ", "")
    a = a.replace("(", "")
    a = a.replace(")", "")
    a = a.replace("**", "*")
    b = b.replace(" ", "")
    b = b.replace("(", "")
    b = b.replace(")", "")
    b = b.replace("**", "*")
    field = a.split("+")
    if len(field) >= 2:
        c = field[1] + "+" + field[0]
    else:
        c = a
    if a == b:
        return True
    elif a.replace("xN", "x2") == b:
        return True
    elif a.replace("xN", "x3") == b:
        return True
    elif a.replace("xN", "x4") == b:
        return True
    elif c.replace("xN", "x2") == b:
        return True
    elif c.replace("xN", "x3") == b:
        return True
    elif c.replace("xN", "x4") == b:
        return True
    elif a.replace("*36+*36", "*36x2") == b:
        return True
    elif a.replace("*68+*68", "*68x2") == b:
        return True
    return False

def main_all_cyp(spec_dir, result_file, cutoff=0):
    gene_class = "CYP"
    truth_dict = load_GeT_RM4()

    spec_result_dict = {} 
    spec_depth_dict = {} 

    
    ## for each folder in the spec_dir, the sample name is the folder name
    for folder in os.listdir(spec_dir):
        ## check if the folder is a folder
        if not os.path.isdir(os.path.join(spec_dir, folder)):
            continue
        sample = folder
        spec_result = os.path.join(spec_dir, folder, f"{folder}.CYP.merge.type.result.txt")
        cyp_2d6_rsult, spec_result_dict[sample], spec_depth_dict[sample] = read_spec_result(spec_result)
        # print (pure_diplotype, spec_result_dict[sample])
    truth_dict = get_shared_sample(truth_dict, spec_result_dict)
    print ("speclong:")

    truth_dict = filter_depth_sample(truth_dict, spec_depth_dict, cutoff)
    spec_result_dict = filter_depth_sample(spec_result_dict, spec_depth_dict, cutoff)

    gene_list, interval_dict =  get_focus_gene(gene_class)
    spec_gene_accuracy_dict = compare_four(truth_dict, spec_result_dict, gene_list, 8, gene_class)
    data = []
    for gene in spec_gene_accuracy_dict:
        data.append(spec_gene_accuracy_dict[gene])

    df = pd.DataFrame(data, columns = ['gene', 'correct', 'total', 'accuracy'])
    df.to_csv(result_file, index=False)

def main_cyp_hprc(pangu_dir, spec_dir, result_file, dataset="1k"):
    if dataset == "1k":
        truth_dict = load_1k_CYP_truth()
    else:
        truth_dict = load_HPRC_CYP_truth()
    pangu_result_dict = {}
    spec_result_dict = {} 
    spec_depth_dict = {} 
    # print (truth_dict)
    # pangu_dir = "/home/wangshuai/00.hla/long/experiments/cyp/cyp_results/pangu_hprc/"
    # spec_dir = "/home/wangshuai/00.hla/long/experiments/cyp/cyp_results/spec_hprc/"

    # pangu_dir = "/home/wangshuai/00.hla/long/experiments/cyp/cyp_results/pangu_hprc_ont/"
    # spec_dir = "/home/wangshuai/00.hla/long/experiments/cyp/cyp_results/spec_hprc_ont/"

    # pangu_dir = "/home/wangshuai/00.hla/long/experiments/cyp/cyp_results/pangu_1k/"
    # spec_dir = "/home/wangshuai/00.hla/long/experiments/cyp/cyp_results/spec_1k/"
    ## for each file with suffix _report.json in the pangu_dir
    for file in os.listdir(pangu_dir):
        if file.endswith("_report.json"):
            sample = file.split("_")[0]
            pangu_result = os.path.join(pangu_dir, file)
            pangu_diplotype = read_pangu_result(pangu_result)
            pangu_result_dict[sample] = pangu_diplotype
    
    #### remove the elements in the truth_dict that are not in the pangu_result_dict
    truth_dict = get_shared_sample(truth_dict, pangu_result_dict)
    print ("pangu", len(pangu_result_dict), "truth", len(truth_dict))
    compare_four(truth_dict, pangu_result_dict, ['CYP2D6'], 8, "CYP")
    
    ## for each folder in the spec_dir, the sample name is the folder name
    for folder in os.listdir(spec_dir):
        ## check if the folder is a folder
        if not os.path.isdir(os.path.join(spec_dir, folder)):
            continue
        sample = folder
        spec_result = os.path.join(spec_dir, folder, f"{folder}.CYP.merge.type.result.txt")
        spec_result_dict[sample], all_gene_result, spec_depth_dict[sample] = read_spec_result(spec_result)
        # print (pure_diplotype, spec_result_dict[sample])
    print ("speclong:")
    compare_four(truth_dict, spec_result_dict, ['CYP2D6'], 8, "CYP")

    cyp_depth_cutoff(truth_dict, spec_depth_dict, spec_result_dict, pangu_result_dict, result_file)



def cyp_depth_cutoff(truth_dict, spec_depth_dict, spec_result_dict, pangu_result_dict, result_file):
    data = []
    cutoff_set = [0, 10, 20, 30, 40, 50]
    #cutoff_set = [x * 5 for x in range(5)]
    #cutoff_set = [x * 10 for x in range(6)]
    for cutoff in cutoff_set:
        print ("###", cutoff)
        cutoff_truth_dict = filter_depth_sample(truth_dict, spec_depth_dict, cutoff)
        cutoff_spec_dict = filter_depth_sample(spec_result_dict, spec_depth_dict, cutoff)
        spec_gene_accuracy_dict = compare_four(cutoff_truth_dict, cutoff_spec_dict, ['CYP2D6'], 8, "CYP")
        data.append(['SpecLong', cutoff] + spec_gene_accuracy_dict['CYP2D6'])

    for cutoff in cutoff_set:
        print ("###", cutoff)
        cutoff_truth_dict = filter_depth_sample(truth_dict, spec_depth_dict, cutoff)
        cutoff_pangu_dict = filter_depth_sample(pangu_result_dict, spec_depth_dict, cutoff)
        pangu_gene_accuracy_dict = compare_four(cutoff_truth_dict, cutoff_pangu_dict, ['CYP2D6'], 8, "CYP")
        data.append(['pangu', cutoff] + pangu_gene_accuracy_dict['CYP2D6'])
    
    # transfrom data to df
    df = pd.DataFrame(data, columns = ['method', 'cutoff', 'gene', 'correct', 'total', 'accuracy'])
    # save to csv
    # df.to_csv("cyp_results/cyp_depth_cutoff.csv", index=False)
    df.to_csv(result_file, index=False)

# def parse_1000g_truth(file):
#     # Region	Population	Sample ID	HLA-A 1	HLA-A 2	HLA-B 1	HLA-B 2	HLA-C 1	HLA-C 2	HLA-DQB1 1	HLA-DQB1 2	HLA-DRB1 1	HLA-DRB1 2
#     # AFR	ACB	HG01879	23:01	68:02	13:02	42:01	08:04	17:01	02:02	04:02	03:02	09:01
#     # AFR	ACB	HG01880	33:03	68:02	40:06	42:01	12:02	17:01	02:01	03:04	03:01	11:06
#     truth_dict = defaultdict(dict)
#     # parse filr use dataframe, but some alleles is format as 03:01/04/09/19, so need to split it to 03:01, 03:04, 03:09, 03:19
#     df = pd.read_csv(file, sep="\t")
#     for index, row in df.iterrows():
#         sample = row['Sample ID']
#         for gene in ['HLA-A', 'HLA-B', 'HLA-C', 'HLA-DQB1', 'HLA-DRB1']:
#             alleles = row[gene].split("/")
#             if len(alleles) == 1:
#                 alleles.append(alleles[0])
#             truth_dict[sample][gene] = [alleles]
    
def parse_1000g_truth(file):
    # Initialize the dictionary to store the results
    truth_dict = defaultdict(dict)
    
    # Read the file into a DataFrame
    df = pd.read_csv(file, sep="\t")
    
    # Process each row in the DataFrame
    for index, row in df.iterrows():
        sample = row['Sample ID']
        
        # Process each gene to extract alleles
        for gene in ['HLA-A 1', 'HLA-A 2', 'HLA-B 1', 'HLA-B 2', 'HLA-C 1', 'HLA-C 2', 'HLA-DQB1 1', 'HLA-DQB1 2', 'HLA-DRB1 1', 'HLA-DRB1 2']:
            gene_idx = int(gene.split(" ")[1]) - 1
            gene_name = gene.split(" ")[0].split("-")[-1]
            
            if gene_name not in truth_dict[sample]:
                truth_dict[sample][gene_name] = [[], []]
            if pd.notna(row[gene]):
                row[gene]=row[gene].strip("*")
                if "/" in row[gene]:
                    alleles = row[gene].split("/")
                    first_field = alleles[0].split(":")[0]
                    for allele in alleles:
                        if first_field in allele:
                            truth_dict[sample][gene_name][gene_idx].append(allele)
                        else:
                            allele_truth = first_field + ":" + allele
                            truth_dict[sample][gene_name][gene_idx].append(allele_truth)
                else:
                    truth_dict[sample][gene_name][gene_idx].append(row[gene])
            else:
                truth_dict[sample][gene_name][gene_idx].append(None)
    
    return truth_dict



def main_1kg_ont_HLA(gene_list, truth_dict, result_dir, gene_class="HLA", step = 2, samples = []):
    ## remove HLA- prefix in gene_list
    if gene_class == "HLA":
        gene_list = [x.split("-")[-1] for x in gene_list]
   
    speclong_result_dict = parse_all_spleclong_1kg_ont_input(gene_class, step, result_dir, samples)
    hla_hla_la_result = parse_all_hlala_input(gene_class, result_dir, samples)
    # print (all_hla_la_result)
    # return
    # all_old_hlala_result = parse_hlala_pacbio()

    # print (new_truth_dict.keys(), new_truth_dict["HG00514.1"].keys())
    # print (all_hla_la_result.keys())
    compare_four(truth_dict, speclong_result_dict, gene_list, 8, gene_class)
    print ("------------------")
    # if gene_class == "IG_TR":
    #     compare_four(new_truth_dict, all_hla_la_result, IG_list, 8, gene_class)
    #     print ("------------------")
    #     compare_four(new_truth_dict, all_hla_la_result, TR_list, 8, gene_class)
    #     print ("------------------")
    
    # result_file = f"{benchmark_result_dir}/{gene_class}/pacbio_{gene_class}.csv"
    # store_results(new_truth_dict, all_hla_la_result, gene_list, result_file)
    # count_report_allele(all_hla_la_result)

def parse_all_spleclong_1kg_ont_input(gene_class, step, outdir, samples):
    all_speclong_result = {}
    if step == 1:
        suffix = ".type.result.txt"
    else:
        suffix = ".final.type.result.txt"
    for sample in samples:
        sample_result = os.path.join(outdir, f"{sample}/{sample}/{sample}.{gene_class}{suffix}")
        if not os.path.exists(sample_result):
            print (sample_result, "does not exist")
            continue
        input_dict = parse_hla_hla_input(sample_result)
        all_speclong_result[sample] = input_dict

def parse_all_hlala_input(gene_class, outdir, samples):
    all_hlala_result = {}
    for sample in samples:
        sample_result = os.path.join(outdir, f"{sample}/{sample}/{sample}/hla/R1_bestguess.txt")
        if not os.path.exists(sample_result):
            print (sample_result, "does not exist")
            continue
        input_dict = parse_hlala_single(sample_result)
        all_hlala_result[sample] = input_dict
    return all_hlala_result

def parse_hlala_single(file_path):
    pass


def read_samples(sample_file):
    samples = []
    with open(sample_file, 'r') as f:
        for line in f:
            samples.append(line.strip())
    return samples


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='compare results')
    # parser.add_argument('truth', help='Input VCF file path')
    # parser.add_argument('input', help='Output VCF file path')
    # parser.add_argument('output', help='Output VCF file path')

    # args = parser.parse_args()
    benchmark_result_dir = "/mnt/d/HLAPro_backup/Nanopore_optimize/benchmark/"

    #### HLA typing in the nanopore amplicon data, comparing with HLA*LA and SpecHLA
    # main()



    # assess_sim()



    #### evaluation in HGSCV2
    # gene_class = "IG_TR"
    # truth_dir = "/mnt/d/HLAPro_backup/Nanopore_optimize/hgscv2_truth_bwa/"
    # result_dir = "/mnt/d/HLAPro_backup/Nanopore_optimize/vdj_results/"

    gene_class = "HLA"
    truth_dir = "/mnt/d/HLAPro_backup/Nanopore_optimize/hgscv2_truth_bwa/"
    result_dir = "/mnt/d/HLAPro_backup/Nanopore_optimize/pacbio_hla/"

    # gene_class = "CYP"
    # truth_dir = "/mnt/d/HLAPro_backup/Nanopore_optimize/hgscv2_truth_bwa/"
    # result_dir = "/mnt/d/HLAPro_backup/Nanopore_optimize/cyp_results/"

    # truth_dir = "/scratch/project/cs_shuaicli/wxd/hla_pacbio_new/hifi/hgscv2_truth_bwa_zip/"
    # result_dir = "/scratch/project/cs_shuaicli/wxd/hla_pacbio_new/hifi/kir_typing_out/"

    truth_dict_1000g=parse_1000g_truth("20181129_HLA_types_full_1000_Genomes_Project_panel.txt")
    print(truth_dict_1000g)


    # step = 2   ### 1 or 2, assess result in step 1 or step 2
    # db_dir = f"../db/{gene_class}/"
    # gene_list, interval_dict =  get_focus_gene(gene_class)
    # gene_mean_len, allele_length_dict = cal_gene_len(db_dir)
    # main_pacbio(gene_list, truth_dir, result_dir, gene_class, step)

    # 1KGP ont HLA
    samples=read_samples("sample.txt")
    step = 2   ### 1 or 2, assess result in step 1 or step 2
    # db_dir = f"../db/{gene_class}/"
    db_dir=f"/gpfs1/scratch/ResearchGroups/cs_shuaicli/wxd/app/SpecLong/db/{gene_class}/"
    gene_list, interval_dict =  get_focus_gene(gene_class)
    gene_mean_len, allele_length_dict = cal_gene_len(db_dir)
    main_1kg_ont_HLA(gene_list, truth_dict_1000g, result_dir, gene_class, step, samples)
    # main_pacbio(gene_list, truth_dir, result_dir, gene_class, step)



    #### TCR evaluation in 11 samples with given truth
    # main_TCR("/mnt/d/HLAPro_backup/Nanopore_optimize/vdj_results_tcr/")
    
    
