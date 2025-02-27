"""
Choose best-mapped allele as reference

wangshuai, wshuai294@gmail.com
"""


import os
import numpy as np
import pickle
import sys
import argparse
import random
import re

# import pulp
# from pulp import LpProblem, LpMinimize, LpMaximize, LpVariable, lpSum, PULP_CBC_CMD, value
from collections import defaultdict
import pysam

from read_objects import My_read
from handle_allele_pair import My_allele_pair
from determine_gene import get_focus_gene, get_folder_list
from db_objects import My_db
from get_allele_depth import Get_depth
from read_binning import filter_fq
from alignment_modules import Read_Type, map2db_blast, map2db
from check_if_homo import if_homo, if_homo2
from downsample_bam import downsample_func
from folder_objects import My_folder

# gene_list = ['A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRB1']



def get_del_CYP2D6():
    CNV_CYP2D6 = ['CYP2D6*68', 'CYP2D6*61', 'CYP2D6*63', 'CYP2D6*4.013', 'CYP2D6*36', 'CYP2D6*83', 'CYP2D6*10', 'CYP2D6*17', 'CYP2D6*13', 'CYP2D6*79','CYP2D6*80', 'CYP2D6*78', 'CYP2D6*67', 'CYP2D6*66', 'CYP2D6*76','CYP2D6*5']
    # RARE_CYP2D6 = ["CYP2D6*65","CYP2D6*147","CYP2D6*101","CYP2D6*146","CYP2D6*171", "CYP2D6*141","CYP2D6*83","CYP2D6*34", "CYP2D6*122", "CYP2D6*121","CYP2D6*86"]
    RARE_CYP2D6 = []
    ignore_CYP2D6 = CNV_CYP2D6 + RARE_CYP2D6
    return  ignore_CYP2D6

def construct_matrix(args, gene, bam, record_candidate_alleles):

    bamfile = pysam.AlignmentFile(bam, 'r')  

    record_read_allele_dict = defaultdict(dict)
    allele_name_dict = defaultdict(int)



    for read in bamfile:
        if read.is_unmapped:
            continue

        my_read = My_read()
        my_read.load_bam(read)

        read_name = my_read.read_name
        allele_name = my_read.allele_name
        gene = my_read.loci_name

        if allele_name not in record_candidate_alleles[gene]:
            continue 

        if allele_name not in record_read_allele_dict[read_name]:
            record_read_allele_dict[read_name][allele_name] = my_read
            allele_name_dict[allele_name] += 1
        # elif record_read_allele_dict[read_name][allele_name].identity < my_read.identity:
        else:
            if args['i'] != 'CYP':
                if record_read_allele_dict[read_name][allele_name].match_num < my_read.match_num:
                    record_read_allele_dict[read_name][allele_name] = my_read
            else:
                if record_read_allele_dict[read_name][allele_name].alignment_score < my_read.alignment_score:
                    record_read_allele_dict[read_name][allele_name] = my_read
            # else:
            #     if my_read.match_num/record_read_allele_dict[read_name][allele_name].match_num > 0.95 and my_read.identity > record_read_allele_dict[read_name][allele_name].identity:
            #         record_read_allele_dict[read_name][allele_name] = my_read
            #     elif record_read_allele_dict[read_name][allele_name].match_num < my_read.match_num:
            #         record_read_allele_dict[read_name][allele_name] = my_read
                # if record_read_allele_dict[read_name][allele_name].identity < my_read.identity:
                #     record_read_allele_dict[read_name][allele_name] = my_read
        # if read_name == "95da644d-6c77-4fab-b63b-5427f480fb95":
        #     if allele_name in ['CYP2D6*1.001', 'CYP2D6*15.003']:
        #         print (allele_name, my_read.match_num, my_read.identity, my_read.alignment_score)
                



    bamfile.close()
    print_read_matrix(args, gene, record_read_allele_dict, allele_name_dict)
    return record_read_allele_dict, allele_name_dict


def construct_matrix_blast(args, gene, blast_file):

    f = open(blast_file, 'r')

    record_read_allele_dict = defaultdict(dict)
    allele_name_dict = defaultdict(int)

    for line in f:
        # skip is empty or comment line
        if not line or line.startswith("#"):
            continue
            
        read_name = line.split("\t")[0]
        allele_name = line.split("\t")[1]
        gene = allele_name.split("*")[0]

        if allele_name not in record_read_allele_dict[read_name]:
            my_read = My_read()
            my_read.load_blast(line)
            record_read_allele_dict[read_name][allele_name] = my_read
            allele_name_dict[allele_name] += 1
        # elif record_read_allele_dict[read_name][allele_name].identity < my_read.identity:
        # elif record_read_allele_dict[read_name][allele_name].match_num < my_read.match_num:
        #     record_read_allele_dict[read_name][allele_name] = my_read
        else:
            record_read_allele_dict[read_name][allele_name].load_second_blast(line) 
    f.close()


    print_read_matrix(args, gene, record_read_allele_dict, allele_name_dict)
    return record_read_allele_dict, allele_name_dict


def print_read_matrix(args, gene, record_read_allele_dict, allele_name_dict):
    allele_name_list = list(allele_name_dict.keys())
    outdir = args["o"] + "/" + args["n"]
    out = open(f"""{my_folder.genes_dir}/{gene}.read.matrix.csv""", 'w')
    first_line = "allele,"
    for allele in allele_name_list:
        first_line += allele + ","
    # first_line += "\n"
    print (first_line, file = out)
    for read_name in record_read_allele_dict:
        line = read_name + ","
        for allele_name in allele_name_list:
            # if read_name == "61babda0-c61c-530c-bdca-b66990689eea" and allele_name == "HLA-B*15:01:01:01":
            #     print (read_name, allele_name, str(record_read_allele_dict[read_name][allele_name].match_num), record_read_allele_dict[read_name][allele_name].identity)
            if allele_name not in record_read_allele_dict[read_name]:
                line += "0/0"+ ","
            else:
                line += str(record_read_allele_dict[read_name][allele_name].match_num) + "/" + str(round(record_read_allele_dict[read_name][allele_name].identity, 4)) + ","
        # line += "\n"
        print (line, file = out)
    out.close()


def model3(gene, record_read_allele_dict, allele_name_dict, record_allele_length):
    ## enumeration
    read_name_list = list(record_read_allele_dict.keys())
    read_num = len(read_name_list)
    allele_name_list = list(allele_name_dict.keys())
    allele_num = len(allele_name_list)


    # print ("read_num:", read_num, "allele_num:", allele_num)

    record_allele_pair_match_len = {}
    record_allele_pair_mismatch = defaultdict(int)
    record_allele_pair_identity = defaultdict(float)
    record_allele_pair_sep_match = {}

    for i in range(allele_num):
        for j in range(i+1, allele_num):

            if args["i"] == "CYP":
                if allele_name_list[i].split(".")[0] in ignore_CYP2D6:
                    continue
                if allele_name_list[j].split(".")[0] in ignore_CYP2D6:
                    continue

            allele_pair_obj = My_allele_pair(allele_name_list[i], allele_name_list[j])
            allele_pair_obj.assign_reads(record_read_allele_dict)

            tag = allele_pair_obj.tag


            allele_pair_obj.allele_1_obj.get_depth(record_allele_length[allele_pair_obj.allele_1])
            allele_pair_obj.allele_2_obj.get_depth(record_allele_length[allele_pair_obj.allele_2])

            allele_pair_obj.allele_1_obj.get_coverage(record_allele_length[allele_pair_obj.allele_1])
            allele_pair_obj.allele_2_obj.get_coverage(record_allele_length[allele_pair_obj.allele_2])

            # if gene in ['HLA-DRB1'] :  # not DPB1 , "HLA-A" "HLA-DPA1", 
            #     depth_cutoff = 0.25
            #     depth_l = [allele_pair_obj.allele_1_obj.depth, allele_pair_obj.allele_2_obj.depth]
            #     if min(depth_l)/max(depth_l) < depth_cutoff:
            #         continue
            
            record_allele_pair_match_len[tag] = allele_pair_obj.pair_obj.match_num
            record_allele_pair_identity[tag] = allele_pair_obj.pair_obj.identity
            # record_allele_pair_identity[tag] = allele_pair_obj.pair_obj.median_identity

            record_allele_pair_sep_match[tag] = {}
            record_allele_pair_sep_match[tag][allele_name_list[i]] = {}
            record_allele_pair_sep_match[tag][allele_name_list[i]]["identity"] = allele_pair_obj.allele_1_obj.identity
            record_allele_pair_sep_match[tag][allele_name_list[i]]["depth"] = allele_pair_obj.allele_1_obj.depth
            record_allele_pair_sep_match[tag][allele_name_list[i]]["coverage"] = allele_pair_obj.allele_1_obj.coverage
            record_allele_pair_sep_match[tag][allele_name_list[j]] = {}
            record_allele_pair_sep_match[tag][allele_name_list[j]]["identity"] = allele_pair_obj.allele_2_obj.identity
            record_allele_pair_sep_match[tag][allele_name_list[j]]["depth"] = allele_pair_obj.allele_2_obj.depth
            record_allele_pair_sep_match[tag][allele_name_list[j]]["coverage"] = allele_pair_obj.allele_2_obj.coverage

    # print ("record_allele_pair_match_len", len(record_allele_pair_match_len))
    tag_list, highest_score = choose_best_alleles(gene, record_allele_pair_match_len, record_allele_pair_identity,record_allele_pair_sep_match)
    first_pair = tag_list[0]
    # print ("first_pair", first_pair)
    tag_list = order_result_pair(tag_list, record_allele_pair_sep_match)
    # print ("tag_list", tag_list)
    type_allele_result =  generate_output(tag_list)          

    return highest_score, type_allele_result, first_pair, record_allele_pair_sep_match

def split_assign_reads(gene, first_pair, record_read_allele_dict, raw_fq):
    allele_name_list = first_pair.split("&")
    # print (allele_name_list)
    allele_pair_obj = My_allele_pair(allele_name_list[0], allele_name_list[1])
    read_assign_dict = allele_pair_obj.assign_reads(record_read_allele_dict)
    # print (read_assign_dict)
    outfile = my_folder.reads_dir + '/%s.fq'%(allele_name_list[0])

    filter_fq(allele_name_list[0], read_assign_dict, raw_fq, outfile)

    outfile = my_folder.reads_dir + '/%s.fq'%(allele_name_list[1])

    filter_fq(allele_name_list[1], read_assign_dict, raw_fq, outfile)

def determine_largest(a, b):
    if a > b:
        return 0
    elif a < b:
        return 1
    else:
        return random.randint(0, 1)

def print_match_results(sorted_record_allele_pair_match_len, record_allele_pair_sep_match, gene, record_allele_pair_identity):
    outdir = args["o"] + "/" + args["n"]
    out = open(f"""{my_folder.genes_dir}/{gene}.allele.match.csv""", 'w')
    print ("index,total_match","total_identity,allele_1,allele_1_identity,allele_1_depth,allele_1_coverage,allele_2,allele_2_identity,allele_2_depth,allele_2_coverage", file = out)

    for i in range(len(sorted_record_allele_pair_match_len)):
        tag = sorted_record_allele_pair_match_len[i][0]
        allele_list = tag.split("&")

        print(
            i, 
            int(sorted_record_allele_pair_match_len[i][1]), 
            round(record_allele_pair_identity[tag],6), 
            allele_list[0],
            round(record_allele_pair_sep_match[tag][allele_list[0]]["identity"],6),
            round(record_allele_pair_sep_match[tag][allele_list[0]]["depth"]),
            round(record_allele_pair_sep_match[tag][allele_list[0]]["coverage"],3),
            allele_list[1], 
            round(record_allele_pair_sep_match[tag][allele_list[1]]["identity"],6), 
            round(record_allele_pair_sep_match[tag][allele_list[1]]["depth"]),
            round(record_allele_pair_sep_match[tag][allele_list[1]]["coverage"],3), 
            sep=",", 
            file=out
        )

    out.close()

def choose_best_alleles(gene, record_allele_pair_match_len, record_allele_pair_identity,record_allele_pair_sep_match):
    # print (record_allele_pair_match_len)
    sorted_record_allele_pair_match_len = sorted(record_allele_pair_match_len.items(), key=lambda x: x[1], reverse=True)
    
    print_match_results(sorted_record_allele_pair_match_len, record_allele_pair_sep_match, gene, record_allele_pair_identity)

    highest_match_score = sorted_record_allele_pair_match_len[0][1]


    # len_diff_cutoff = 1e-2  #1e-4 
    len_diff_cutoff = 1e-3 
    ide_diff_cutoff = 1e-4   # 1e-4  # 0.001
    if gene in ["HLA-B"]:
        len_diff_cutoff = 1e-4

    # if gene  in ["HLA-DPA1"]:
    #     len_diff_cutoff = 1e-2
    if gene in ["HLA-DRB1"]:
        len_diff_cutoff = 0.1 ## 5e-2
        ide_diff_cutoff = 2e-4
    # if gene  in ["HLA-C"]:
    #     len_diff_cutoff =  1e-3
    if args['i'] == "CYP":
        len_diff_cutoff = 0.1
        ide_diff_cutoff = 1e-5
        
    if args['i'] == "KIR":
        len_diff_cutoff = 0.05
        ide_diff_cutoff = 1e-5
    if gene == 'KIR3DP1':
        len_diff_cutoff = 0.1
    if gene == 'KIR3DL3':
        len_diff_cutoff = 0.01
    if gene == 'KIR2DP1':
        len_diff_cutoff = 0.001
    if gene == 'KIR2DS4':
        len_diff_cutoff = 0.005

    # if gene  in ["DPB1"]:
    #     len_diff_cutoff = 1e-2

    good_length_dict = {}
    for i in range(len(sorted_record_allele_pair_match_len)):
        tag = sorted_record_allele_pair_match_len[i][0]
        allele_list = tag.split("&")
        # print (allele_list)

        # print(
        #     "#pass len", int(sorted_record_allele_pair_match_len[i][1]),
        #     round(record_allele_pair_identity[tag],3),
        #     allele_list[0],
        #     round(record_allele_pair_sep_match[tag][allele_list[0]]["identity"],3),
        #     round(record_allele_pair_sep_match[tag][allele_list[0]]["depth"]),
        #     round(record_allele_pair_sep_match[tag][allele_list[0]]["coverage"],3),
        #     allele_list[1],
        #     round(record_allele_pair_sep_match[tag][allele_list[1]]["identity"],3),
        #     round(record_allele_pair_sep_match[tag][allele_list[1]]["depth"]),
        #     round(record_allele_pair_sep_match[tag][allele_list[1]]["coverage"],3),
        #     sep="\t"
        # )
        if args['i'] == "KIR":
            if abs(record_allele_pair_sep_match[tag][allele_list[0]]["identity"] - record_allele_pair_sep_match[tag][allele_list[1]]["identity"]) > 0.1:
                continue
            depth_list = [record_allele_pair_sep_match[tag][allele_list[0]]["depth"], record_allele_pair_sep_match[tag][allele_list[1]]["depth"]]
            if max(depth_list) != 0 and min(depth_list)/max(depth_list) < 0.2:
                continue

        if (highest_match_score - sorted_record_allele_pair_match_len[i][1])/highest_match_score <= len_diff_cutoff:
            good_length_dict[tag] = record_allele_pair_identity[tag]
        else:
            break
    print ("# allele pairs that pass match len cutoff", len(good_length_dict))

    if len(good_length_dict) == 0:
        print ("# no allele pair selected, use less strict criteria", len(good_length_dict))
        for i in range(len(sorted_record_allele_pair_match_len)):
            tag = sorted_record_allele_pair_match_len[i][0]
            allele_list = tag.split("&")
            if (highest_match_score - sorted_record_allele_pair_match_len[i][1])/highest_match_score <= len_diff_cutoff:
                good_length_dict[tag] = record_allele_pair_identity[tag]
            else:
                break

    identity_sorted_list = sorted(good_length_dict.items(), key=lambda x: x[1], reverse=True)
    match_len_with_max_identity = identity_sorted_list[0][1]
    full_result_list = []
    for i in range(len(identity_sorted_list)):
        if (match_len_with_max_identity - identity_sorted_list[i][1])/match_len_with_max_identity <= ide_diff_cutoff:
            print ("identity", identity_sorted_list[i])
            full_result_list.append(identity_sorted_list[i][0])
            
            if len(full_result_list) >= 30:
                break

    return full_result_list, highest_match_score

def get_most_common_allele(allele_list, allele):
    max_digit_num = 0
    for a in allele_list:
        same_digit_num = cal_sim_of_alleles(allele, a)
        if same_digit_num > max_digit_num:
            max_digit_num = same_digit_num
    return max_digit_num

def generate_output(tag_list):
    type_allele_result = [[], []]
    for tag in tag_list:
        pair = tag.split("&")
        if len(type_allele_result[0]) == 0:
            for i in range(2):
                if pair[i] not in type_allele_result[i]:
                    type_allele_result[i].append(pair[i])
        else:
            max_digit_num_list1 = [get_most_common_allele(type_allele_result[0], pair[0]), get_most_common_allele(type_allele_result[1], pair[0])]
            max_digit_num_list2 = [get_most_common_allele(type_allele_result[1], pair[1]), get_most_common_allele(type_allele_result[1], pair[1])]
            # print (max_digit_num_list1, max_digit_num_list2)

            if max_digit_num_list1[0] > max_digit_num_list1[1]:
                direction = "forward"
            elif max_digit_num_list1[0] < max_digit_num_list1[1]:
                direction = "backward"
            elif max_digit_num_list2[0] > max_digit_num_list2[1]:
                direction = "backward"
            elif max_digit_num_list2[0] < max_digit_num_list2[1]:
                direction = "forward"
            else:
                direction = "forward"
            # print (direction)
            if direction == "forward":
                for i in range(2):
                    if pair[i] not in type_allele_result[i]:
                        type_allele_result[i].append(pair[i])
            else:
                for i in range(2):
                    if pair[1-i] not in type_allele_result[i]:
                        type_allele_result[i].append(pair[1-i])
        
        # print (pair, type_allele_result)
    for i in range(2):
        type_allele_result[i] = ",".join(list(type_allele_result[i]))
    # print ("type_allele_result", type_allele_result)
    return type_allele_result

def order_result_pair(type_allele_result, record_allele_pair_sep_match):
    sort_type_allele_result = []
    for tag in type_allele_result:
        allele_list = tag.split("&")

        if record_allele_pair_sep_match[tag][allele_list[1]]["depth"] > record_allele_pair_sep_match[tag][allele_list[0]]["depth"]:
            tag = allele_list[1] + "&" + allele_list[0]

        # if record_allele_pair_sep_match[tag][allele_list[1]]["identity"] > record_allele_pair_sep_match[tag][allele_list[0]]["identity"]:
        #     tag = allele_list[1] + "&" + allele_list[0]

        sort_type_allele_result.append(tag)
    # print ("sort_type_allele_result", sort_type_allele_result)
    return sort_type_allele_result

def cal_sim_of_alleles(allele1, allele2):
    clean_allele1 = allele1.split("*")[1].split(":")      
    clean_allele2 = allele2.split("*")[1].split(":")   
    same_digit_num = 0
    for i in range(min([len(clean_allele1), len(clean_allele2)])): 
        if clean_allele1[i] == clean_allele2[i]:
            same_digit_num += 2
        else:
            break
    return same_digit_num

def output_spechla_format(args, result_dict):
    outdir = args["o"] + "/" + args["n"]
    out = open(f"{outdir}/hla.new.result.txt", 'w')
    print ("#", file = out)
    print ("sample", end = "\t", file = out)
    for gene in gene_list:
        for i in range(2):
            print (f"{gene}_{i+1}", end = "\t", file = out)

    print (f"\n{args['n']}", end = "\t", file = out)
    for gene in gene_list:
        if len(result_dict[gene]) == 1:
            result_dict[gene].append(result_dict[gene][0])
 
        for i in range(2):
            print (result_dict[gene][i], end = "\t", file = out)

    # print (f"\nMatch_length", end = "\t", file = out)
    # for gene in gene_list:
    #     if len(result_dict[gene]) == 1:
    #         result_dict[gene].append(result_dict[gene][0])
    #     for i in range(2):
    #         print (allele_match_dict[result_dict[gene][i]], end = "\t", file = out)

    # print (f"\nObjective value:\t", end = "\t", file = out)
    out.close()

def output_hlala_format(args, result_dict, reads_num_dict, homo_p_value_dict, p_value_cutoff=0.0001):
    outdir = args["o"] + "/" + args["n"]
    result = f"""{outdir}/{args["n"]}.{args["i"]}.type.result.txt"""
    f = open(result, 'w')
    print (my_db.version_info, file = f)
    print ("Locus   Chromosome      Allele  Reads_num   homo_flag   Homo_p  Hete_pair", file = f)
    for gene in gene_list:
        # if len(result_dict[gene]) == 1:
        #     result_dict[gene].append(result_dict[gene][0])
        homo_flag = False
        if homo_p_value_dict[gene] != 'NA' and homo_p_value_dict[gene] < p_value_cutoff:
            homo_flag = True
        
        if not homo_flag:
            print (gene, "hete", result_dict[gene][0], "****", result_dict[gene][1], "\n\n")
        else:
            print (gene, "homo", result_dict[gene][0], "****", result_dict[gene][1], "\n\n")

        for ch in [1, 2]:
            if re.search('HLA-C\*04:01:01:11', result_dict[gene][ch-1]):
                result_dict[gene][ch-1] += ",HLA-C*04:01:01:01"
            result_dict[gene][ch-1] = result_dict[gene][ch-1].replace(',', ';')
            type_allele = result_dict[gene][ch-1]
            if homo_flag:  # homo
                type_allele = result_dict[gene][0]
            print (gene, ch, type_allele, reads_num_dict[gene], homo_flag, homo_p_value_dict[gene], result_dict[gene][ch-1], sep="\t", file = f)
    f.close()
    print ("result is", result)

def select_candidate_allele(record_allele_length, allele_match_dict, allele_mismatch_dict, max_allele_num):
    ## define dict to save the raw depth of each allele
    allele_raw_dp = {}
    allele_raw_identity = {}

    for allele in allele_match_dict:
        allele_raw_dp[allele] = allele_match_dict[allele]/record_allele_length[allele]
        allele_raw_identity[allele] = allele_match_dict[allele]/(allele_match_dict[allele] + allele_mismatch_dict[allele])
    filter_allele_name_dict = {}
    ## sort allele_raw_dp based on value
    sorted_dict = sorted(allele_raw_dp.items(), key=lambda x: x[1], reverse=True)
    # sorted_dict = sorted(allele_raw_identity.items(), key=lambda x: x[1], reverse=True)
    for i in range(min([max_allele_num, len(sorted_dict)])):
        print (sorted_dict[i], allele_match_dict[sorted_dict[i][0]])
        filter_allele_name_dict[sorted_dict[i][0]] = 1
    return filter_allele_name_dict
    

def cal_allele_match_len(record_read_allele_dict):
    allele_match_dict = defaultdict(int)
    allele_mismatch_dict = defaultdict(int)
    for read_name in record_read_allele_dict:
        for allele_name in record_read_allele_dict[read_name]:
            allele_match_dict[allele_name] += record_read_allele_dict[read_name][allele_name].match_num
            allele_mismatch_dict[allele_name] += record_read_allele_dict[read_name][allele_name].mismatch_num
    return allele_match_dict, allele_mismatch_dict

def main():
    result_dict = {}
    reads_num_dict = {}
    homo_p_value_dict = {}


    for gene in gene_list:
        print (f"start type {gene} for {args['n']}...\n")

        gene_fq=f"{my_folder.reads_dir}/{gene}.long_read.fq.gz"
        my_db.get_allele_length(gene)
        record_allele_length = my_db.allele_len_dict


        #  load alignment from bam
        bam, depth_file, sort_depth_file = map2db(args, gene, my_db, my_folder, args["max_read_num"],args["align_method"])
        get_depth = Get_depth(depth_file)
        get_depth.record_depth()
        record_candidate_alleles, record_allele_length_no_use = get_depth.select(sort_depth_file, gene_list, args["candidate_allele_num"])
        # record_candidate_alleles[gene].add("HLA-C*06:02:01:01")

        if args["align_method"] == "minimap2":
            record_read_allele_dict, allele_name_dict = construct_matrix(args, gene, bam, record_candidate_alleles)

            ### select the alleles with top match num
            # allele_match_dict, allele_mismatch_dict = cal_allele_match_len(record_read_allele_dict)
            # allele_name_dict = select_candidate_allele(record_allele_length, allele_match_dict, allele_mismatch_dict, 100)

        elif args["align_method"] == "blastn":
            #  load alignment from blast
            blast_file = map2db_blast(args, gene, my_db)
            record_read_allele_dict, allele_name_dict = construct_matrix_blast(args, gene, blast_file)





            allele_name_dict = {}
            for allele in record_candidate_alleles[gene]:
                allele_name_dict[allele] = 1

        reads_num_dict[gene] = len(record_read_allele_dict)
        print ("finish matrix construction")
        


        
        print (gene, "read num:", len(record_read_allele_dict), "mapped allele num:", len(allele_name_dict))
        if len(record_read_allele_dict) >= args["min_read_num"] and len(allele_name_dict) > 1:
            objective_value, type_allele_result, first_pair, record_allele_pair_sep_match = model3( gene, record_read_allele_dict, allele_name_dict, record_allele_length)
            split_assign_reads(gene, first_pair, record_read_allele_dict, gene_fq)  # output assigned reads to fastq

            print (gene, type_allele_result, "\n\n")
            homo_p_value = if_homo2(record_allele_pair_sep_match, first_pair)
            homo_p_value_dict[gene] = homo_p_value
            # if homo_p_value < 0.001:
            #     type_allele_result = [type_allele_result[0]] 


            

            # homo_hete_ratio = (objective_value - allele_match_dict[type_allele_result[0].split(",")[0]])/ allele_match_dict[type_allele_result[0].split(",")[0]]
            
            # homo_hete_ratio_cutoff = args["b"]

            # if gene == "HLA-DPA1":
            #     homo_hete_ratio_cutoff = 0.001
            # # if gene == "DQA1":
            # #     homo_hete_ratio_cutoff = 0.01

            # print (homo_hete_ratio, homo_hete_ratio_cutoff)

            # if gene == "HLA-C":
            #     if cal_sim_of_alleles(type_allele_result[0].split(",")[0], type_allele_result[1].split(",")[0]) != 6:
            #         if homo_hete_ratio  <  homo_hete_ratio_cutoff:
            #             type_allele_result = [type_allele_result[0]]
            
            # else:
            #     if homo_hete_ratio  <  homo_hete_ratio_cutoff:
            #         type_allele_result = [type_allele_result[0]]      
        else:
            print ("support read number is too low or the mapped allele num is too low for %s, skip typing"%(gene))
            type_allele_result  = ['-', '-']
            homo_p_value_dict[gene] = 'NA'


        result_dict[gene] = type_allele_result

    # output_spechla_format(args, result_dict)
    output_hlala_format(args, result_dict, reads_num_dict, homo_p_value_dict, args["hete_p"])
    


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="HLA Typing with long-read data.", add_help=False, \
    usage="python3 %(prog)s -h", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    required = parser.add_argument_group("Required arguments")
    optional = parser.add_argument_group("Optional arguments")
    # required.add_argument("-f", type=str, help="IMGT reference.", metavar="\b")
    required.add_argument("-r", type=str, help="Long-read fastq file. PacBio or Nanopore.", metavar="\b")
    required.add_argument("-n", type=str, help="Sample ID", metavar="\b")
    required.add_argument("-o", type=str, help="The output folder to store the typing results.", metavar="\b", default="./output")
    required.add_argument("-i", type=str, help="HLA,KIR,CYP",metavar="\b", default="HLA")
    optional.add_argument("-j", type=int, help="Number of threads.", metavar="\b", default=5)
    optional.add_argument("--candidate_allele_num", type=int, help="Maintain this number of alleles for ILP step.", metavar="\b", default=100)
    optional.add_argument("--min_read_num", type=int, help="min support read number for each locus.", metavar="\b", default=2)
    optional.add_argument("--max_read_num", type=int, help="max support read number for each locus.", metavar="\b", default=1000)
    optional.add_argument("-b", type=float, help="The match length increase ratio lower than this value is homo [0-1].", metavar="\b", default=0.0007)
    optional.add_argument("--hete_p", type=float, help="Hete pvalue.", metavar="\b", default=1e-30)
    optional.add_argument("--db", type=str, help="db dir.", metavar="\b", default=sys.path[0] + "/../db/")
    # optional.add_argument("-g", type=int, help="Whether use G group resolution annotation [0|1].", metavar="\b", default=0)
    # optional.add_argument("-m", type=int, help="1 represents typing, 0 means only read assignment", metavar="\b", default=1)
    optional.add_argument("-y", type=str, help="Read type, [nanopore|pacbio|pacbio-hifi].", metavar="\b", default="pacbio")
    optional.add_argument("--align_method", type=str, help="[minimap2|blastn|bwa].", metavar="\b", default="minimap2")
    optional.add_argument("--test", type=bool, help="for test.", metavar="\b", default=False)
    # optional.add_argument("--max_depth", type=int, help="maximum depth for each HLA locus. Downsample if exceed this value.", metavar="\b", default=10000)
    # optional.add_argument("-u", type=str, help="Choose full-length or exon typing. 0 indicates full-length, 1 means exon.", metavar="\b", default="0")
    optional.add_argument("-rt", "--RNA_type", type=str, help="traditional,2D,Direct,SIRV",metavar="\b", default="traditional")
    optional.add_argument("--seq_tech", type=str, help="Amplicon sequencing or WGS sequencing [wgs|amplicon|rna].", metavar="\b", default="wgs")
    optional.add_argument("-h", "--help", action="help")
    args = vars(parser.parse_args()) 

    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(0)


    ignore_CYP2D6 = get_del_CYP2D6()
    # gene_list, interval_dict =  get_focus_gene(args)
    my_db = My_db(args)
    my_folder = My_folder(args)
    # db_folder=os.path.dirname(my_db.full_cds_db) if args["seq_tech"] == "rna" else os.path.dirname(my_db.full_db)
    db_folder = os.path.dirname(my_db.full_db)
    gene_list = get_folder_list(db_folder)
    if args["test"]:
        # gene_list = ['HLA-A', 'HLA-B', 'HLA-C', 'HLA-DPA1', 'HLA-DPB1', 'HLA-DQA1', 'HLA-DQB1', 'HLA-DRB1']
        gene_list = ['HLA-DRB1']
    main()
