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

# import pulp
# from pulp import LpProblem, LpMinimize, LpMaximize, LpVariable, lpSum, PULP_CBC_CMD, value
from collections import defaultdict
import pysam

from read_objects import My_read
from handle_allele_pair import My_allele_pair
from filter_reads import examine_reads
from determine_gene import get_focus_gene
from db_objects import My_db
from get_allele_depth import Get_depth
from read_binning import filter_fq
# gene_list = ['A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRB1']
  


def construct_matrix(args, gene, bam, record_candidate_alleles, record_allele_length):

    bamfile = pysam.AlignmentFile(bam, 'r')  

    record_read_allele_dict = defaultdict(dict)
    allele_name_dict = defaultdict(int)



    for read in bamfile:
        if read.is_unmapped:
            continue

        my_read = My_read(read)

        read_name = read.query_name
        allele_name = read.reference_name
        gene = my_read.loci_name

        if allele_name not in record_candidate_alleles[gene]:
            continue 

        if allele_name not in record_read_allele_dict[read_name]:
            record_read_allele_dict[read_name][allele_name] = my_read
            allele_name_dict[allele_name] += 1
        # elif record_read_allele_dict[read_name][allele_name].identity < my_read.identity:
        elif record_read_allele_dict[read_name][allele_name].match_num < my_read.match_num:
            record_read_allele_dict[read_name][allele_name] = my_read



    bamfile.close()
    print_read_matrix(args, gene, record_read_allele_dict, allele_name_dict)
    return record_read_allele_dict, allele_name_dict


def print_read_matrix(args, gene, record_read_allele_dict, allele_name_dict):
    allele_name_list = list(allele_name_dict.keys())
    outdir = args["o"] + "/" + args["n"]
    out = open(f"""{outdir}/{args["n"]}.{gene}.read.matrix.csv""", 'w')
    first_line = "allele,"
    for allele in allele_name_list:
        first_line += allele + ","
    # first_line += "\n"
    print (first_line, file = out)
    for read_name in record_read_allele_dict:
        line = read_name + ","
        for allele_name in allele_name_list:
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
            allele_pair_obj = My_allele_pair(allele_name_list[i], allele_name_list[j])
            allele_pair_obj.assign_reads(record_read_allele_dict)

            tag = allele_pair_obj.tag


            allele_pair_obj.allele_1_obj.get_depth(record_allele_length[allele_pair_obj.allele_1])
            allele_pair_obj.allele_2_obj.get_depth(record_allele_length[allele_pair_obj.allele_2])

            if gene in ["HLA-DPA1", 'HLA-DRB1'] :  # not DPB1
                depth_cutoff = 0.25
                depth_l = [allele_pair_obj.allele_1_obj.depth, allele_pair_obj.allele_2_obj.depth]
                if min(depth_l)/max(depth_l) < depth_cutoff:
                    continue
            
            record_allele_pair_match_len[tag] = allele_pair_obj.pair_obj.match_num
            record_allele_pair_identity[tag] = allele_pair_obj.pair_obj.identity

            record_allele_pair_sep_match[tag] = {}
            record_allele_pair_sep_match[tag][allele_name_list[i]] = {}
            record_allele_pair_sep_match[tag][allele_name_list[i]]["identity"] = allele_pair_obj.allele_1_obj.identity
            record_allele_pair_sep_match[tag][allele_name_list[i]]["depth"] = allele_pair_obj.allele_1_obj.depth
            record_allele_pair_sep_match[tag][allele_name_list[j]] = {}
            record_allele_pair_sep_match[tag][allele_name_list[j]]["identity"] = allele_pair_obj.allele_2_obj.identity
            record_allele_pair_sep_match[tag][allele_name_list[j]]["depth"] = allele_pair_obj.allele_2_obj.depth

    # print ("record_allele_pair_match_len", len(record_allele_pair_match_len))
    tag_list, highest_score = choose_best_alleles(gene, record_allele_pair_match_len, record_allele_pair_identity,record_allele_pair_sep_match)
    tag_list = order_result_pair(tag_list, record_allele_pair_sep_match)
    type_allele_result =  generate_output(tag_list)          



    return highest_score, type_allele_result, tag_list[0]

def split_assign_reads(gene, first_pair, record_read_allele_dict, outdir, raw_fq):
    allele_name_list = first_pair.split("&")
    # print (allele_name_list)
    allele_pair_obj = My_allele_pair(allele_name_list[0], allele_name_list[1])
    read_assign_dict = allele_pair_obj.assign_reads(record_read_allele_dict)
    # print (read_assign_dict)
    outfile = outdir + '/%s.fq'%(allele_name_list[0])

    filter_fq(allele_name_list[0], read_assign_dict, raw_fq, outfile)

    outfile = outdir + '/%s.fq'%(allele_name_list[1])

    filter_fq(allele_name_list[1], read_assign_dict, raw_fq, outfile)


def determine_largest(a, b):
    if a > b:
        return 0
    elif a < b:
        return 1
    else:
        return random.randint(0, 1)

def select_alleles_by_identity_diff(record_allele_pair_sep_match, identity_cutoff=0.07, depth_cutoff= 0.25):
    allele_pair_after_identity_diff = {}
    for tag in record_allele_pair_sep_match:
        allele_list = tag.split("&")
        if abs(record_allele_pair_sep_match[tag][allele_list[1]]["identity"]  - record_allele_pair_sep_match[tag][allele_list[0]]["identity"]) > identity_cutoff:
            continue
        depth_l = [record_allele_pair_sep_match[tag][allele_list[0]]["depth"], record_allele_pair_sep_match[tag][allele_list[1]]["depth"]]
        if min(depth_l)/max(depth_l) < depth_cutoff:
            continue
        allele_pair_after_identity_diff[tag] = 1
    print ("alleles_after_identity_diff pair", len(allele_pair_after_identity_diff))
    return allele_pair_after_identity_diff

def print_match_results(sorted_record_allele_pair_match_len, record_allele_pair_sep_match, gene, record_allele_pair_identity):
    outdir = args["o"] + "/" + args["n"]
    out = open(f"""{outdir}/{args["n"]}.{gene}.allele.match.csv""", 'w')
    print ("index,total_match","total_identity,allele_1,allele_1_identity,allele_1_depth,allele_2,allele_2_identity,allele_2_depth", file = out)

    data = []
    for i in range(len(sorted_record_allele_pair_match_len)):
        tag = sorted_record_allele_pair_match_len[i][0]
        allele_list = tag.split("&")
        print (i, int(sorted_record_allele_pair_match_len[i][1]), round(record_allele_pair_identity[tag],3), allele_list[0], \
            round(record_allele_pair_sep_match[tag][allele_list[0]]["identity"],3),round(record_allele_pair_sep_match[tag][allele_list[0]]["depth"]),\
             allele_list[1], round(record_allele_pair_sep_match[tag][allele_list[1]]["identity"],3), round(record_allele_pair_sep_match[tag][allele_list[1]]["depth"]), sep = ",", file = out)
    out.close()

def choose_best_alleles(gene, record_allele_pair_match_len, record_allele_pair_identity,record_allele_pair_sep_match):
    # print (record_allele_pair_match_len)
    sorted_record_allele_pair_match_len = sorted(record_allele_pair_match_len.items(), key=lambda x: x[1], reverse=True)
    
    print_match_results(sorted_record_allele_pair_match_len, record_allele_pair_sep_match, gene, record_allele_pair_identity)

    highest_match_score = sorted_record_allele_pair_match_len[0][1]


    len_diff_cutoff = 1e-4 
    ide_diff_cutoff = 0.001
    # if gene  in ["DQA1", "DRB1", "DPA1"]:
    #     len_diff_cutoff = 1e-3
    if gene  in ["HLA-DPA1"]:
        len_diff_cutoff = 1e-2
    if gene  in ["HLA-DRB1"]:
        len_diff_cutoff = 5e-2
    # if gene  in ["DPB1"]:
    #     len_diff_cutoff = 1e-2

    good_length_dict = {}
    for i in range(len(sorted_record_allele_pair_match_len)):
        tag = sorted_record_allele_pair_match_len[i][0]
        allele_list = tag.split("&")
        print ("#pass len", int(sorted_record_allele_pair_match_len[i][1]), round(record_allele_pair_identity[tag],3), allele_list[0], \
            round(record_allele_pair_sep_match[tag][allele_list[0]]["identity"],3),round(record_allele_pair_sep_match[tag][allele_list[0]]["depth"]),\
             allele_list[1], round(record_allele_pair_sep_match[tag][allele_list[1]]["identity"],3), round(record_allele_pair_sep_match[tag][allele_list[1]]["depth"]), sep="\t")
        if (highest_match_score - sorted_record_allele_pair_match_len[i][1])/highest_match_score <= len_diff_cutoff:
            good_length_dict[tag] = record_allele_pair_identity[tag]
        else:
            break

    # for i in range(len(sorted_record_allele_pair_match_len)):
    #     tag = sorted_record_allele_pair_match_len[i][0]
    #     allele_list = tag.split("&")
    #     print ("#", int(sorted_record_allele_pair_match_len[i][1]), round(record_allele_pair_identity[tag],3), allele_list[0], \
    #         round(record_allele_pair_sep_match[tag][allele_list[0]]["identity"],3),round(record_allele_pair_sep_match[tag][allele_list[0]]["depth"]),\
    #          allele_list[1], round(record_allele_pair_sep_match[tag][allele_list[1]]["identity"],3), round(record_allele_pair_sep_match[tag][allele_list[1]]["depth"]), sep="\t")

    #     if i > 100:
    #         break
    identity_sorted_list = sorted(good_length_dict.items(), key=lambda x: x[1], reverse=True)
    match_len_with_max_identity = identity_sorted_list[0][1]
    full_result_list = []
    for i in range(len(identity_sorted_list)):
        if (match_len_with_max_identity - identity_sorted_list[i][1])/match_len_with_max_identity <= ide_diff_cutoff:
            print ("identity", identity_sorted_list[i])
            full_result_list.append(identity_sorted_list[i][0])
            
            if len(full_result_list) >= 30:
                break


    # first consider identity, then consider match length
    # if gene == "DQA1":
    #     len_diff_cutoff = 1e-2
    #     ide_diff_cutoff = 4e-3
    #     identity_sorted_list = sorted(record_allele_pair_identity.items(), key=lambda x: x[1], reverse=True)
    #     match_len_with_max_identity = identity_sorted_list[0][1]
    #     good_identity_dict = {}
    #     full_result_list = []
    #     for i in range(len(identity_sorted_list)):
    #         if (match_len_with_max_identity - identity_sorted_list[i][1])/match_len_with_max_identity <= ide_diff_cutoff:
    #             good_identity_dict[identity_sorted_list[i][0]] = record_allele_pair_match_len[identity_sorted_list[i][0]]
    #             print (identity_sorted_list[i], record_allele_pair_match_len[identity_sorted_list[i][0]])
    #         else:
    #             break
    #     match_len_sorted_list = sorted(good_identity_dict.items(), key=lambda x: x[1], reverse=True)
    #     match_len_with_max_identity = match_len_sorted_list[0][1]
    #     full_result_list = []
    #     for i in range(len(match_len_sorted_list)):
    #         if (match_len_with_max_identity - match_len_sorted_list[i][1])/match_len_with_max_identity <= len_diff_cutoff:
    #             print ("match length", match_len_sorted_list[i])
    #             full_result_list.append(match_len_sorted_list[i][0])
                
    #         if len(full_result_list) >= 30:
    #             break



    # print ("\nhomo\n")
    # for i in range(len(sorted_record_allele_pair_match_len)):
    #     tag = sorted_record_allele_pair_match_len[i][0]
    #     if tag.split("&")[0] == tag.split("&")[1]:
    #         print (tag.split("&"), int(sorted_record_allele_pair_match_len[i][1]), record_allele_pair_identity[tag], sep="\t")

    # print ("\ntruth\n")
    # tag = 'DQA1*01:03:01:01'  + "&" + "DQA1*01:04:02"
    # print (tag, record_allele_pair_match_len[tag], record_allele_pair_identity[tag], sep="\t")

    

    return full_result_list, highest_match_score


def generate_output(tag_list):
    type_allele_result = ['', '']
    for tag in tag_list:
        pair = tag.split("&")
        if len(type_allele_result[0]) == 0:
            type_allele_result = pair
            if type_allele_result[0] == 'C*04:01:01:11':
                type_allele_result[0] = 'C*04:01:01:01,C*04:01:01:11'
            if type_allele_result[1] == 'C*04:01:01:11':
                type_allele_result[1] = 'C*04:01:01:01,C*04:01:01:11'
        else:
            
            for p in pair:
                max_digit_num = [0, 0]
                keep = True
                for i in range(2):
                    for allele in type_allele_result[i].split(","):
                        if p == allele: 
                            keep = False
                        same_digit_num = cal_sim_of_alleles(p, allele)
                        # print (p, allele, same_digit_num)
                        if same_digit_num > max_digit_num[i]:
                            max_digit_num[i] = same_digit_num
                if keep:
                    if p == 'C*04:01:01:11':
                        p = 'C*04:01:01:01,C*04:01:01:11'
                    if max_digit_num[0] > max_digit_num[1]:
                        type_allele_result[0] += "," + p
                    elif max_digit_num[0] < max_digit_num[1]:
                        type_allele_result[1] += "," + p
                    else:
                        type_allele_result[0] += "," + p
                        type_allele_result[1] += "," + p

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

def map2db(args, gene):

    minimap_para = ''
    if args["y"] == "pacbio":
        minimap_para = " -x map-pb "
    elif args["y"] == "nanopore":
        minimap_para = " -x map-ont "

    outdir = args["o"] + "/" + args["n"]
    sam = outdir + "/" + args["n"] + "." + gene + ".db.sam"
    bam = outdir + "/" + args["n"] + "." + gene + ".db.bam"
    depth_file = outdir + "/" + args["n"] + "." + gene + ".db.depth"
    sort_depth_file = outdir + "/" + args["n"] + "." + gene + ".db.sort.depth.txt"
    # ref={args["f"] }
    # map raw reads to database

    ref = my_db.get_gene_all_alleles(gene)
    # ref="/mnt/d/HLAPro_backup/Nanopore_optimize/SpecHLA/db/HLA/whole/HLA_A.fasta"

    alignDB_order = f"""
    fq={args["r"]}
    
    outdir={args["o"]}/{args["n"]}
    sample={args["n"]}

    fq=$outdir/{gene}.long_read.fq.gz

    ref={ref}

    minimap2 -t {args["j"]} {minimap_para} -p 0.1 -N 100000 -a $ref $fq > {sam}

    # bwa index $ref
    # bwa mem -R '@RG\\tID:foo\\tSM:bar' -a -t {args["j"]} $ref $fq > {sam}

    samtools view -bS -F 0x800  {sam} | samtools sort - >{bam}
    samtools index {bam}
    samtools depth -aa {bam}>{depth_file}
    rm {sam}
    echo alignment done.
    """
    # if the depth_file is not detected 
    if not os.path.exists(depth_file):
        os.system(alignDB_order)
    else:
        print("Depth file is detected.")

    return bam, depth_file, sort_depth_file

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

def output_hlala_format(args, result_dict, reads_num_dict):
    outdir = args["o"] + "/" + args["n"]
    result = f"""{outdir}/{args["n"]}.{args["i"]}.type.result.txt"""
    f = open(result, 'w')
    # print ("#", version_info, file = f)
    print ("Locus   Chromosome      Allele  Reads_num", file = f)
    for gene in gene_list:
        if len(result_dict[gene]) == 1:
            result_dict[gene].append(result_dict[gene][0])
        for ch in [1, 2]:
            result_dict[gene][ch-1] = result_dict[gene][ch-1].replace(',', ';')
            print (gene, ch, result_dict[gene][ch-1], reads_num_dict[gene], sep="\t", file = f)
    f.close()
    print ("result is", result)

def main(args):

    if not os.path.exists(args["o"]):
        os.system("mkdir %s"%(args["o"]))
    outdir = args["o"] + "/" + args["n"]
    if not os.path.exists(outdir):
        os.system("mkdir %s"%(outdir))
    # outdir = args["o"]

    result_dict = {}
    reads_num_dict = {}


    for gene in gene_list:
        print (f"start type {gene} for {args['n']}...\n")
        allele_match_dict = defaultdict(int)
        allele_read_num_dict = defaultdict(int)
        gene_fq=f"{outdir}/{gene}.long_read.fq.gz"

        bam, depth_file, sort_depth_file = map2db(args, gene)

        get_depth = Get_depth(depth_file)
        get_depth.record_depth()
        record_candidate_alleles, record_allele_length = get_depth.select(sort_depth_file, gene_list, args["candidate_allele_num"])
        # gene="A"
        # print ("record_candidate_alleles num:", len(record_candidate_alleles[gene]))
        print (bam)

        record_read_allele_dict, allele_name_dict = construct_matrix(args, gene, bam, record_candidate_alleles, record_allele_length)
        print ("finish matrix construction")
        reads_num_dict[gene] = len(record_read_allele_dict)
        # record_read_allele_dict = examine_reads(record_read_allele_dict)
        # print (len(record_read_allele_dict)) 

        for read_name in record_read_allele_dict:
            for allele_name in record_read_allele_dict[read_name]:
                allele_match_dict[allele_name] += record_read_allele_dict[read_name][allele_name].match_num
                allele_read_num_dict[allele_name] += 1
        
        # for allele_name in allele_match_dict:
        #     print ("speclong depth", allele_name, allele_match_dict[allele_name], round(allele_match_dict[allele_name]/record_allele_length[allele_name]),allele_read_num_dict[allele_name])

        print (gene, "read num:", len(record_read_allele_dict), "mapped allele num:", len(allele_name_dict))
        if len(record_read_allele_dict) >= args["min_read_num"] and len(allele_name_dict) > 1:
            objective_value, type_allele_result, first_pair = model3( gene, record_read_allele_dict, allele_name_dict, record_allele_length)
            split_assign_reads(gene, first_pair, record_read_allele_dict, outdir, gene_fq)  # output assigned reads to fastq

            print (gene, type_allele_result, "\n\n")
            

            homo_hete_ratio = (objective_value - allele_match_dict[type_allele_result[0].split(",")[0]])/ allele_match_dict[type_allele_result[0].split(",")[0]]
            
            homo_hete_ratio_cutoff = args["b"]

            if gene == "HLA-DPA1":
                homo_hete_ratio_cutoff = 0.001
            # if gene == "DQA1":
            #     homo_hete_ratio_cutoff = 0.01

            print (homo_hete_ratio, homo_hete_ratio_cutoff)

            if gene == "HLA-C":
                if cal_sim_of_alleles(type_allele_result[0].split(",")[0], type_allele_result[1].split(",")[0]) != 6:
                    if homo_hete_ratio  <  homo_hete_ratio_cutoff:
                        type_allele_result = [type_allele_result[0]]
            
            else:
                if homo_hete_ratio  <  homo_hete_ratio_cutoff:
                    type_allele_result = [type_allele_result[0]]      
        else:
            print ("support read number is too low or the mapped allele num is too low for %s, skip typing"%(gene))
            type_allele_result  = ['-', '-']


        result_dict[gene] = type_allele_result
        print (gene, type_allele_result, "\n\n")
    # output_spechla_format(args, result_dict)
    output_hlala_format(args, result_dict, reads_num_dict)
    


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
    optional.add_argument("--min_read_num", type=int, help="min support read number.", metavar="\b", default=2)
    optional.add_argument("-b", type=float, help="The match length increase ratio lower than this value is homo [0-1].", metavar="\b", default=0.0007)
    optional.add_argument("--db", type=str, help="db dir.", metavar="\b", default=sys.path[0] + "/../db/")
    # optional.add_argument("-g", type=int, help="Whether use G group resolution annotation [0|1].", metavar="\b", default=0)
    # optional.add_argument("-m", type=int, help="1 represents typing, 0 means only read assignment", metavar="\b", default=1)
    optional.add_argument("-y", type=str, help="Read type, [nanopore|pacbio].", metavar="\b", default="pacbio")
    # optional.add_argument("-u", type=str, help="Choose full-length or exon typing. 0 indicates full-length, 1 means exon.", metavar="\b", default="0")
    optional.add_argument("-h", "--help", action="help")
    args = vars(parser.parse_args()) 

    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(0)



    gene_list, interval_dict =  get_focus_gene(args)
    my_db = My_db(args)

    # gene_list = ["HLA-A"]
    main(args)
