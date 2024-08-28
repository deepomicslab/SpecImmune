"""
For each HGT gene sequence and the blast result
assign the suitable types for it

It is perfect choice if the max-match-len allele and the max-identity allele are same.
Otherwise, we need to balance the match length and identity to chose the suitable alleles.

wangshuai, Apr 23, 2024
"""

import os
import re
import sys
import pysam
import argparse
from collections import defaultdict

from determine_gene import get_focus_gene, get_folder_list
from db_objects import My_db
from folder_objects import My_folder
from select_best_reference_alleleV2 import get_del_CYP2D6


def get_1_element(lst):
    return lst[1]

def get_2_element(lst):
    return lst[2]

def get_3_element(lst):
    return lst[3]

def change_allele_name(raw, new):
    with open(raw, "r") as infile, open(new, "w") as outfile:
        for line in infile:
            if line.startswith(">"):
                header = line.strip()[1:]
                contig_name = header.split()[1]
                new_header = f">{contig_name}\n"
                outfile.write(new_header)
            else:
                outfile.write(line)

def resort_list_with_same_alleles(sorted_list, first_index, second_index):
    flag = True
    while flag:
        flag = False
        new_sorted_list = sorted_list.copy()
        for i in range(len(sorted_list) - 1):
            if sorted_list[i][first_index] == sorted_list[i+1][first_index] and sorted_list[i+1][second_index] > sorted_list[i][second_index]:
                new_sorted_list[i] = sorted_list[i+1]
                new_sorted_list[i+1] = sorted_list[i]
                flag = True
        sorted_list = new_sorted_list.copy()
    # print (sorted_list[:5])
    return sorted_list
    
def get_max_alleles(sorted_list, index):
    max_value = sorted_list[0][index]
    max_allele_list = []
    for list in sorted_list:
        if list[index] == max_value:
            # max_allele_list.append(list[0])
            list = [str(x) for x in list]
            max_allele_list.append(">".join(list))
        else:
            break
    return max_allele_list

def extract_four_digits(full_name):
    a = full_name.split("*")[1]
    array = a.split(":")
    return array[0] + ":" + array[1]

def compare_match_len_and_identity(match_sorted_list, identity_sorted_list, truth_alleles):
    max_match_len = match_sorted_list[0][1]
    match_len_with_max_identity = identity_sorted_list[0][1]

    max_identity = identity_sorted_list[0][3]
    identiy_with_max_match_len = match_sorted_list[0][3]

    match_len_diff_ratio = (max_match_len - match_len_with_max_identity) / match_len_with_max_identity
    identity_diff_ratio = (max_identity - identiy_with_max_match_len) / identiy_with_max_match_len

    print ("match_len_diff_ratio", match_len_diff_ratio, "identity_diff_ratio", identity_diff_ratio)
    get_help_from_1000G = False

    if extract_four_digits(match_sorted_list[0][0]) in truth_alleles and extract_four_digits(identity_sorted_list[0][0]) not in truth_alleles:
        select_allele_list = match_sorted_list[0]
        get_help_from_1000G = True
    elif extract_four_digits(match_sorted_list[0][0]) not in truth_alleles and extract_four_digits(identity_sorted_list[0][0]) in truth_alleles:
        select_allele_list = identity_sorted_list[0]
        get_help_from_1000G = True
    elif identiy_with_max_match_len < 0.999:
        select_allele_list = identity_sorted_list[0]
    elif match_len_diff_ratio < identity_diff_ratio:
        select_allele_list = identity_sorted_list[0]
    elif match_len_diff_ratio < 0.3:
        select_allele_list = identity_sorted_list[0]
    # elif identity_diff_ratio < 0.005:
    #     select_allele_list = match_sorted_list[0]
    else:
        print (" no determine")
        
    # if get_help_from_1000G == False:
    print ("check to determine use highest identity or match length in person.")
    for allele_info in match_sorted_list[:5]:
        print(allele_info)
    print ("match bases**************************")

    
    for allele_info in identity_sorted_list[:5]:
        print(allele_info)
    print ("identity **************************")
    for allele_info in identity_sorted_list:
        if allele_info[0] == "DRB1*16:02:01:03":
            print (allele_info)
    
    print ("selected allele is ", select_allele_list[0])
    return select_allele_list
    

def select_by_alignment(align_list, gene):
    if len(align_list) == 0:
        return [], []
    full_result_list = []
    best_match_result = []
    # match_sorted_list = sorted(align_list, key=get_1_element, reverse = True)
    # if gene == "HLA-A":
    #     print("align:", align_list[0])
    match_sorted_list = sorted(align_list, key=get_2_element, reverse = True)  # match length - mismatch
    match_sorted_list = resort_list_with_same_alleles(match_sorted_list, 1, 3)
    identity_sorted_list = sorted(align_list, key=get_3_element, reverse = True)
    identity_sorted_list = resort_list_with_same_alleles(identity_sorted_list, 3, 1)
    max_match_len_alleles = get_max_alleles(match_sorted_list, 2)
    max_identity_alleles = get_max_alleles(identity_sorted_list, 3)

    # f = open(all_align_result_file, 'w')
    # for mat in match_sorted_list:
    #     print (mat, file = f)
    # f.close()
    len_diff_cutoff = 0.02
    ide_diff_cutoff = 0.0004

    # if gene == 'KIR3DP1':
    #     len_diff_cutoff = 0.15
    # print ("ienti:",identity_sorted_list)
    intersection_alleles = list(set(max_match_len_alleles) & set(max_identity_alleles))   
    # print (">>>>>>>>>", match_sorted_list[:10])
    if len(intersection_alleles) > 0:
        for z in range(len(intersection_alleles)):
            select_allele_list = intersection_alleles[z].split(">")
            full_result_list.append(select_allele_list)  
        best_match_result = full_result_list.copy() 
        if gene == "HLA-DPB1" or gene == "HLA-C":  
            full_result_list = []
            for i in range(len(identity_sorted_list)):
                if (identity_sorted_list[0][3] - identity_sorted_list[i][3])/identity_sorted_list[0][3] <= len_diff_cutoff:
                    full_result_list.append(identity_sorted_list[i])
                    if len(full_result_list) >= 10:
                        break
    else:
        # 
        # print (">>>>>>>>>>max_match_len allele and max_identity allele don't match, report possible alleles")  
        max_match_len = match_sorted_list[0][2]
        match_len_with_max_identity = identity_sorted_list[0][1]

        good_length_list = []
        for i in range(len(match_sorted_list)):
            if (max_match_len - match_sorted_list[i][2])/max_match_len <= len_diff_cutoff:
                good_length_list.append(match_sorted_list[i])
                # if gene=="HLA-A":
                #     print("match_sorted_list", match_sorted_list[i])
        # print (len(good_length_list), "len(good_length_list)")
        identity_sorted_list = sorted(good_length_list, key=get_3_element, reverse = True)
        identity_sorted_list = resort_list_with_same_alleles(identity_sorted_list, 3, 1)
        match_len_with_max_identity = identity_sorted_list[0][3]
        # full_result_list = identity_sorted_list[:20]
        for i in range(len(identity_sorted_list)):
            if (match_len_with_max_identity - identity_sorted_list[i][3])/match_len_with_max_identity <= ide_diff_cutoff:
                full_result_list.append(identity_sorted_list[i])
                # print (match_len_with_max_identity, identity_sorted_list[i][3], (match_len_with_max_identity - identity_sorted_list[i][3])/match_len_with_max_identity)
                if len(full_result_list) >= 15:
                    break

    return full_result_list, best_match_result

def get_read_num_from_step1(step1_result):
    # check if the step1_result file exists
    if not os.path.exists(step1_result):
        print(f"{step1_result} does not exist")
        return 0
    read_num_dict = {}
    with open(step1_result, "r") as f:
        for line in f:
            if line.startswith("Locus"):
                continue
            if line.startswith("#"):
                continue
            field = line.strip().split()
            locus = field[0]
            hap = field[1]
            allele = field[2]
            tag = f"{locus}_{hap}"
            read_num = int(field[3])
            read_num_dict[tag] = [read_num, allele]
    return read_num_dict


def output(record_best_match, record_one_allele, gene_list, result_file, version_info, record_all_match,read_num_dict):
    f = open(result_file, 'w')
    print (version_info, file = f)
    print ("Locus\tChromosome\tGenotype\tMatch_info\tReads_num\tStep1_type\tOne_guess", file = f)
    for gene in gene_list:
        for ch in [1, 2]:
            tag = f"{gene}_{ch}"
            alleles = ''
            alleles_info = ''
            read_num, step1_alleles = read_num_dict[tag]
            for a in record_best_match[gene][ch]:
                # print (a)
                alleles += a[0] + ";"
                info = record_all_match[tag][a[0]]
                alleles_info += f"{a[0]}|{int(info[2])}|{round(info[3],6)};"
            out_alleles = alleles[:-1]
            alleles_info = alleles_info[:-1]
            if out_alleles == "":
                out_alleles = "NA"
                alleles_info = "NA"
            if re.search('HLA-C\*04:01:01:11', out_alleles):
                out_alleles += ";HLA-C*04:01:01:01"
            print (gene, ch, out_alleles,alleles_info, read_num,step1_alleles, record_one_allele[gene][ch], sep="\t", file = f)
    f.close()

def get_blast_info(blastfile, tag, record_all_match):
    # print("blastfile", blastfile)
    ## check if blastfile exists
    if not os.path.exists(blastfile):
        print(f"{blastfile} does not exist")
        sys.exit(0)
    align_list = []
    for line in open(blastfile):
        field = line.strip().split()
        if field[0] != tag:
            continue
        if args["i"] == "CYP":
            if field[1].split(".")[0] in ignore_CYP2D6:
                continue
        align_info = [field[1], int(field[3]), int(field[3])-float(field[2]), 1-float(field[2])/int(field[3]), 'x', 0, 0]
        record_all_match[tag][field[1]] = align_info
        # if field[1] in ["HLA-A*03:08:01:02", "HLA-A*03:01:01:12"]:
        #     print("field", field, align_info)    
        align_list.append(align_info)
    # print("align_list_0", align_list)
    identity_sorted_list = sorted(align_list, key=get_3_element, reverse = True)
    return identity_sorted_list, record_all_match

def get_digit_array(allele):
    array = allele.split("*")[1].split(":")
    for i in range(len(array)):
        array[i] = array[i].lstrip('0')
        if array[i][-1] in ['N', 'L', 'S', 'C', 'A' ,'Q']:
            array[i] = float('inf')
        
        array[i] = float(array[i])
    add_num = 4 - len(array)
    for i in range(add_num):
        array.append(float('inf'))
    # print (array)
    return array

def compare_num(a, b):
    if args["i"] == "CYP":
        num1 = float(a.split("*")[1])
        num2 = float(b.split("*")[1])
        if num1 <= num2:
            return True
        else:
            return False
    else:
        array1 = get_digit_array(a)
        array2 = get_digit_array(b)
        for i in range(4):
            if array1[i] < array2[i]:
                return True
            elif array1[i] > array2[i]:
                return False
            else:
                pass
        return True
                

def sort_allele_by_number(allele_list):
    all_right = False
    while not all_right:
        all_right = True
        for i in range(1, len(allele_list)):
            right_order = compare_num(allele_list[i-1], allele_list[i])
            # print (allele_list[i-1], allele_list[i], right_order)
            all_right = all_right and right_order
            if not right_order:
                allele_list[i-1], allele_list[i] = allele_list[i], allele_list[i-1]
                # break
    # print ("sort", allele_list)
    return allele_list

def report_one_allele(full_result_list, best_match_result):
    consider_result = []
    if len(best_match_result) > 0:
        for a in best_match_result:
            consider_result.append(a[0])
    else:
        for a in full_result_list:
            consider_result.append(a[0])

    # print (consider_result)

    if len(consider_result) == 0:
        return "NA"
    elif len(consider_result) == 1:
        return consider_result[0]
    else:
        if args['i'].strip() == "CYP":# or args['i'].strip() == "KIR":
            return consider_result[0]
        else:
            sort_consider_result = sort_allele_by_number(consider_result)
            return sort_consider_result[0]

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Refine typing.", add_help=False, \
    usage="python3 %(prog)s -h", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    required = parser.add_argument_group("Required arguments")
    optional = parser.add_argument_group("Optional arguments")
    # required.add_argument("-r", type=str, help="Long-read fastq file. PacBio or Nanopore.", metavar="\b")
    # required.add_argument("-1", type=str, help="Assembly file of the first haplotype in fasta formate", metavar="\b")
    # required.add_argument("-2", type=str, help="Assembly file of the second haplotype in fasta formate", metavar="\b")
    required.add_argument("-n", type=str, help="Sample ID", metavar="\b")
    required.add_argument("-o", type=str, help="The output folder to store the typing results.", metavar="\b", default="./output")
    optional.add_argument("--db", type=str, help="db dir.", metavar="\b", default=sys.path[0] + "/../db/")
    required.add_argument("-i", type=str, help="HLA,KIR,CYP",metavar="\b", default="HLA")
    optional.add_argument("-rt", "--RNA_type", type=str, help="traditional,2D,Direct,SIRV",metavar="\b", default="traditional")
    optional.add_argument("--seq_tech", type=str, help="Amplicon sequencing or WGS sequencing [wgs|amplicon|rna].", metavar="\b", default="wgs")

    
    # optional.add_argument("-g", type=int, help="Whether use G group resolution annotation [0|1].", metavar="\b", default=0)
    # optional.add_argument("-u", type=str, help="Choose full-length or exon typing. 0 indicates full-length, 1 means exon.", metavar="\b", default="0")
    optional.add_argument("-h", "--help", action="help")
    args = vars(parser.parse_args()) 

    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(0)

    ignore_CYP2D6 = get_del_CYP2D6()
    my_db = My_db(args)    
    my_folder = My_folder(args)
    # gene_list, interval_dict =  get_focus_gene(args)

    # db_folder=os.path.dirname(my_db.full_cds_db) if args["seq_tech"] == "rna" else os.path.dirname(my_db.full_db)
    db_folder = os.path.dirname(my_db.full_db)
    gene_list = get_folder_list(db_folder)
    
    # result_file = result_path + "/" + "hlala.like.results.txt"
    result_file = f"{my_folder.sample_prefix}.{args['i']}.final.type.result.txt"
    step1_result = f"{my_folder.sample_prefix}.{args['i']}.type.result.txt"

    record_best_match, record_one_allele = {}, {}
    record_all_match = defaultdict(dict) ## for output
    blastfile = f"{my_folder.sequence_dir}/{args['i']}.blast.summary.txt"

    print (f"optimize typing results by balancing alignment length and identity ")
    for hap_index in range(2):
        for gene in gene_list:
            if gene not in record_best_match:
                record_best_match[gene] = {}
                record_one_allele[gene] = {}
            tag = f"{gene}_{hap_index+1}"
            align_list, record_all_match = get_blast_info(blastfile, tag, record_all_match)
            full_result_list, best_match_result = select_by_alignment(align_list, gene)
            record_best_match[gene][hap_index+1] = full_result_list   
            one_allele = report_one_allele(full_result_list, best_match_result)
            # print ("one allele", one_allele)
            record_one_allele[gene][hap_index+1] = one_allele
    
    read_num_dict = get_read_num_from_step1(step1_result)
    output(record_best_match, record_one_allele, gene_list, result_file, my_db.version_info, record_all_match,read_num_dict)
    print (f"The refined typing results is in {result_file}")