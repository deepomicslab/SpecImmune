import os
import numpy as np
import pickle
import sys
import argparse
import random
import re
import gzip
import pysam
import sys
from collections import defaultdict

from db_objects import My_db
from determine_gene import get_focus_gene
from bed_objects import Bed_db

def cal_gene_depth(depth_file, gene_interval_dict):
    gene_depth_dict = {}
    gene_pos_dict = {}
    with open(depth_file, "r") as f:
        for line in f:
            line = line.strip().split("\t")
            chrom = line[0]
            pos = int(line[1])
            depth = int(line[2])
            for gene in gene_interval_dict[chrom]:
                start = gene_interval_dict[chrom][gene][0]
                end = gene_interval_dict[chrom][gene][1]
                gene_pos_dict[gene] = [chrom] + list(gene_interval_dict[chrom][gene])
                if start <= pos <= end:
                    if gene not in gene_depth_dict:
                        gene_depth_dict[gene] = []
                    gene_depth_dict[gene].append(depth)
                    break

    gene_mean_depth_dict = {}
    for gene in gene_depth_dict:
        gene_mean_depth_dict[gene] = round(np.mean(gene_depth_dict[gene]),1)
        # print (gene, gene_mean_depth_dict[gene], gene_pos_dict[gene])

    return gene_mean_depth_dict, gene_pos_dict

def fast_cal_gene_depth(depth_file, gene_interval_dict):
    gene_pos_dict = {}
    store_all_depth = defaultdict(list)
    with open(depth_file, "r") as f:
        for line in f:
            line = line.strip().split("\t")
            chrom = line[0]
            pos = int(line[1])
            depth = int(line[2])
            store_all_depth[chrom].append(depth)

    gene_mean_depth_dict = {}

    for chrom in store_all_depth:
        segment_start = 0
        segment_end = float("inf")
        if chrom not in gene_interval_dict:

            ## for chr:start-end
            # field = chrom.split(":")
            # pure_chrom = field[0]
            # segment_start = int(field[1].split("-")[0])
            # segment_end = int(field[1].split("-")[1])

            field = chrom.split("_")
            pure_chrom = field[0]
            segment_start = int(field[1])
            segment_end = int(field[2])

        else:
            pure_chrom = chrom

        for gene in gene_interval_dict[pure_chrom]:


            #### check if the gene is in the segment
            pure_start = gene_interval_dict[pure_chrom][gene][0]
            pure_end = gene_interval_dict[pure_chrom][gene][1]
            if pure_end < segment_start or pure_start > segment_end:
                continue

            start = pure_start - segment_start
            end = pure_end - segment_start
            gene_pos_dict[gene] = [chrom, start, end] 
            # if gene == "IGHD3-22":
            #     print (gene, chrom, start, end, store_all_depth[chrom][start-1:end])
            gene_mean_depth_dict[gene] = round(np.mean(store_all_depth[chrom][start-1:end]),1)




    return gene_mean_depth_dict, gene_pos_dict


def parse_phase_blocks(vcf_path):
    # Open the VCF file
    vcf_file = pysam.VariantFile(vcf_path)

    phase_blocks = {}
    
    for record in vcf_file.fetch():
        for sample in record.samples:
            sample_data = record.samples[sample]
            if 'PS' in sample_data:
                phase_set = sample_data['PS']
                chrom = record.chrom
                pos = record.pos

                if phase_set not in phase_blocks:
                    phase_blocks[phase_set] = {
                        'chrom': chrom,
                        'start': pos,
                        'end': pos
                    }
                else:
                    phase_blocks[phase_set]['end'] = pos

    return phase_blocks

def parse_gene_phase_blocks(vcf_path):
    # Open the VCF file
    vcf_file = pysam.VariantFile(vcf_path)

    blocks = set()
    variant_num = 0
    hete_variant_num = 0
    
    for record in vcf_file.fetch():
        variant_num += 1
        for sample in record.samples:
            sample_data = record.samples[sample]
            if 'PS' in sample_data:
                phase_set = sample_data['PS']
                chrom = record.chrom
                pos = record.pos
                # print (phase_set, chrom, pos, sample_data)
                if phase_set:
                    blocks.add(str(phase_set))
            if 'UG' in sample_data:
                unphased_genotype = sample_data['UG']
                if unphased_genotype == "0/1" or unphased_genotype == "1/2":
                    hete_variant_num += 1
    blocks = list(blocks)
    # print (blocks, variant_num)
    if len(blocks) > 0:
        my_block = ";".join(blocks)
    else:
        my_block = "NA"
    return my_block, [variant_num,hete_variant_num]

## given phase_blocks, estimate the phase set of each gene
def assign_phase_block(phase_blocks, gene_interval_dict):
    gene_phase_dict = {}
    for gene in gene_interval_dict:
        chrom = gene_interval_dict[gene][0]
        start = gene_interval_dict[gene][1]
        end = gene_interval_dict[gene][2]
        for phase_set, interval in phase_blocks.items():
            if interval['chrom'] == chrom and start <= interval['start'] <= end:
                gene_phase_dict[gene] = phase_set
                break
    # print (gene_phase_dict)
    return gene_phase_dict

def load_raw_result(raw_result, gene_mean_depth_dict, gene_phase_dict, new_result, min_depth,focus_gene_list):
    """
    raw result is like:
    sample  gene    allele  score   len     start   end     chr     hap
    IG_TR_dp50_acc98_2      IGHA1   IGHA1*02        100.000 393     1025372 1025764 chr14_igh       hap1
    IG_TR_dp50_acc98_2      IGHA2   IGHA2*02        100.000 393     905141  905533  chr14_igh       hap1
    IG_TR_dp50_acc98_2      IGHD    IGHD*01 100.000 324     1159103 1159426 chr14_igh       hap1
    IG_TR_dp50_acc98_2      IGHD2-21        IGHD2-21*02     100.000 28      1206752 1206779 chr14_igh       hap1
    IG_TR_dp50_acc98_2      IGHE    IGHE*04 100.000 330     918270  918599  chr14_igh       hap1
    """  
    # construct new result
    store_raw = defaultdict(dict)
    sample = ''
    with open(raw_result, "r") as f:
        # skip the header
        header = f.readline()
        for line in f:
            line = line.strip().split("\t")
            sample = line[0]
            gene = line[1]
            hap = line[-1]
            store_raw[gene][hap] = line
    
    with open(new_result, "w") as f:
        header = ["sample", "gene", "depth", "phase_set", "allele_1", "score_1", "length_1", "hap_1", "chrom_1", "start_1", "end_1",\
                    "allele_2", "score_2", "length_2", "hap_2", "chrom_2", "start_2", "end_2"]
        print (*header, sep="\t", file=f)
        # for gene in store_raw:
        for gene in my_db.order_gene_list:
            if gene not in focus_gene_list:
                continue
            raw_has_content = True
            if gene in gene_mean_depth_dict:
                depth = gene_mean_depth_dict[gene]
            else:
                depth = "NA"
            if gene in gene_phase_dict:
                phase = gene_phase_dict[gene]
            else:
                phase = "NA"
            if "hap1" in store_raw[gene]:
                sample, gene, allele_1, score_1, length_1, start_1, end_1, chrom_1, hap_1 = store_raw[gene]["hap1"]
            else:
                raw_has_content = False
                allele_1, score_1, length_1, start_1, end_1, chrom_1, hap_1 = "NA", "NA", "NA", "NA", "NA", "NA", "NA"
            if "hap2" in store_raw[gene]:
                sample, gene, allele_2, score_2, length_2, start_2, end_2, chrom_2, hap_2 = store_raw[gene]["hap2"]
            else:
                raw_has_content = False
                allele_2, score_2, length_2, start_2, end_2, chrom_2, hap_2 = "NA", "NA", "NA", "NA", "NA", "NA", "NA"
            
            if depth == "NA" or depth < min_depth:
                allele_1, score_1, length_1, hap_1, chrom_1, start_1, end_1 = "NA", "NA", "NA", "NA", "NA", "NA", "NA"
                allele_2, score_2, length_2, hap_2, chrom_2, start_2, end_2 = "NA", "NA", "NA", "NA", "NA", "NA", "NA"
            if raw_has_content:
                print (sample, gene, depth, phase, allele_1, score_1, length_1, hap_1, chrom_1, start_1, end_1,\
                        allele_2, score_2, length_2, hap_2, chrom_2, start_2, end_2, sep="\t", file=f)

            # print (sample, gene, depth, phase, allele_1, score_1, length_1, hap_1, chrom_1, start_1, end_1,\
            #         allele_2, score_2, length_2, hap_2, chrom_2, start_2, end_2, sep="\t")

def read_blast_result(blast_file, len_cutoff = 0.9):
    # allele, identity, length = "NA", "NA", "NA"
    result = ["NA", "NA", "NA"]
    with open(blast_file) as in_file:
        for line in in_file:
            if line.startswith("#"):
                continue
            # Split the line into a list of values
            values = line.split()
            # Assign values to variables, converting numerical values to integers
            allele, identity, length, mis, gap, start, end, score = values[1], float(values[2]), int(values[3]), int(values[4]), int(values[5]), int(values[6]), int(values[7]), float(values[11])

            gene = allele.split("*")[0]
              
            if result[0] == "NA":
                result = [allele, identity, length]
            else:
                # if length >= result[2] - 5 and identity > result[1]:
                if length/result[2] > len_cutoff and identity > result[1]:
                    result = [allele, identity, length]
                if gene == "IGHV4-34" and identity > result[1]:
                    result = [allele, identity, length]

    return result

def get_consensus(phased_vcf, focus_gene_list, tmp_dir, args, depth_bed_file, gene_mean_depth_dict, short_cutoff = 50):
    
    interval_dict = bed_db.get_gene_interval(bed_db.lite_gene_file)
    all_result_dict = {}
    gene_phase_dict = {}
    variant_num_dict = {}
    for gene in focus_gene_list:
        # if gene != "IGHD5-5":
        #     continue
        if gene not in interval_dict:
            print (f"{gene} not in the gene interval dict")
            continue
        if gene not in gene_mean_depth_dict or gene_mean_depth_dict[gene] <= 0:
            # print (f"{gene} not in the gene depth dict")
            continue
        gene_hg38_len = hg38_gene_info[gene][3]

        if gene_hg38_len < short_cutoff:
            para = " -word_size 4 -task blastn-short "
        else:
            para = ''

        all_result_dict[gene] = {"hap1":{}, "hap2":[]}
        for hap in range(1, 3):
            cmd = f"samtools faidx {args['hg38']} {interval_dict[gene]} | bcftools consensus -H {hap} --mask {depth_bed_file} {phased_vcf} > {tmp_dir}/{gene}.{hap}.fasta"
            os.system(cmd)
            blast_db = my_db.get_blast_index(gene)
            cmd = f"blastn {para} -query {tmp_dir}/{gene}.{hap}.fasta -db {blast_db} -outfmt 7 -max_target_seqs 3000 -num_threads {args['j']} > {tmp_dir}/{gene}.{hap}.blast.txt"
            # print (cmd)
            os.system(cmd)
            result = read_blast_result(f"{tmp_dir}/{gene}.{hap}.blast.txt")
            all_result_dict[gene][f"hap{hap}"] = result
        cmd = f"""bcftools view {phased_vcf} -r {interval_dict[gene]} -Oz -o {tmp_dir}/{gene}.phase.vcf.gz
                    tabix -f {tmp_dir}/{gene}.phase.vcf.gz"""
        os.system(cmd)
        my_block, variant_num = parse_gene_phase_blocks(f"{tmp_dir}/{gene}.phase.vcf.gz")
        gene_phase_dict[gene] = my_block
        variant_num_dict[gene] = variant_num
        # break
    return all_result_dict, gene_phase_dict, variant_num_dict

def generate_output_file(all_result_dict, gene_mean_depth_dict, gene_phase_dict, new_result, min_depth, focus_gene_list, variant_num_dict):
    
    with open(new_result, "w") as f:
        header = ["gene", "depth", "phase_set", "allele_1", "score_1", "length_1", \
                  "hap_1", "allele_2", "score_2", "length_2", "hap_2","hg38_chrom", "hg38_len", "variant_num","hete_variant_num"]
        print (*header, sep="\t", file=f)
        # for gene in store_raw:
        for gene in my_db.order_gene_list:
            if gene not in focus_gene_list:
                continue
            if gene not in hg38_gene_info:
                print (f"{gene} not in the hg38_gene_info")
                continue
            if gene not in all_result_dict:
                continue
            raw_has_content = True
            if gene in gene_mean_depth_dict:
                depth = gene_mean_depth_dict[gene]
            else:
                depth = "NA"
            if gene in gene_phase_dict:
                phase = gene_phase_dict[gene]
            else:
                phase = "NA"
            if "hap1" in all_result_dict[gene]:
                allele_1, score_1, length_1 = all_result_dict[gene]["hap1"]
            else:
                allele_1, score_1, length_1 = "NA", "NA", "NA"

            if "hap2" in all_result_dict[gene]:
                allele_2, score_2, length_2 = all_result_dict[gene]["hap2"]
            else:
                allele_2, score_2, length_2 = "NA", "NA", "NA"
            
            if depth == "NA" or depth <= 0:
                allele_1, score_1, length_1, hap_1, chrom_1, start_1, end_1 = "NA", "NA", "NA", "NA", "NA", "NA", "NA"
                allele_2, score_2, length_2, hap_2, chrom_2, start_2, end_2 = "NA", "NA", "NA", "NA", "NA", "NA", "NA"
            else:
                print (gene, depth, phase, allele_1, score_1, length_1, "hap1", allele_2, score_2, length_2, "hap2", \
                        hg38_gene_info[gene][0], hg38_gene_info[gene][3], variant_num_dict[gene][0], variant_num_dict[gene][1], sep="\t", file=f)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate depth and assign phased block to IG_TR genes.", add_help=False, \
    usage="python3 %(prog)s -h", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    required = parser.add_argument_group("Required arguments")
    optional = parser.add_argument_group("Optional arguments")
    # required.add_argument("--phased_vcf", type=str, help="phased_vcf.", metavar="\b")
    # required.add_argument("--depth_file", type=str, help="depth_file", metavar="\b")
    # required.add_argument("--raw_result", type=str, help="raw_result", metavar="\b")
    required.add_argument("-n", type=str, help="Sample ID", metavar="\b")
    required.add_argument("-o", type=str, help="The output folder to store the typing results.", metavar="\b", default="./output")
    required.add_argument("-i", type=str, help="IG_TR",metavar="\b", default="IG_TR")
    optional.add_argument("--db", type=str, help="db dir.", metavar="\b", default=sys.path[0] + "/../db/")
    optional.add_argument("--hg38", type=str, help="referece fasta file, used by IG_TR typing, generated by extract_VDJ_segments_from_hg38.py", metavar="\b", default=sys.path[0] + "/../VDJ_ref/IG_TR.segment.fa")
    optional.add_argument("-k", type=int, help="The mean depth in a window lower than this value will be masked by N, set 0 to avoid masking", metavar="\b", default=5)
    optional.add_argument("-j", type=int, help="Number of threads.", metavar="\b", default=5)
    optional.add_argument("-h", "--help", action="help")
    args = vars(parser.parse_args()) 

    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(0)
    
    my_db = My_db(args)
    bed_db = Bed_db()
    hg38_gene_info = bed_db.get_hg38_gene_interval() 
    focus_gene_list, xx =  get_focus_gene("IG_TR")

    # for chrom in my_db.gene_interval_dict:
    #     for gene in my_db.gene_interval_dict[chrom]:
    #         print (gene, chrom, my_db.gene_interval_dict[chrom][gene])

    depth_file = args["o"] + "/" + args["n"] + "/" + args["n"] + ".depth.txt"
    depth_bed_file = args["o"] + "/" + args["n"] + "/low_depth.bed"
    phased_vcf = args["o"] + "/" + args["n"] + "/" + args["n"] + ".phase.norm.vcf.gz"
    raw_result = args["o"] + "/" + args["n"] + "/" + args["n"] + ".IG.TR.allele.txt"
    new_result = args["o"] + "/" + args["n"] + "/" + args["n"] + ".IG_TR_typing_result.txt"
    tmp_dir = args["o"] + "/" + args["n"] + "/genes/"
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir) 
    # print (my_db.gene_interval_dict)
    gene_mean_depth_dict, gene_pos_dict = fast_cal_gene_depth(depth_file, my_db.gene_interval_dict)
    # get_phased_block(args["phased_vcf"])


    # phase_blocks = parse_phase_blocks(phased_vcf)
    
    # for phase_set, interval in phase_blocks.items():
    #     print(f"Phase Set {phase_set}: Chromosome {interval['chrom']}, Start {interval['start']}, End {interval['end']}")
    # gene_phase_dict = assign_phase_block(phase_blocks, gene_pos_dict)
    all_result_dict, gene_phase_dict, variant_num_dict = get_consensus(phased_vcf, focus_gene_list, tmp_dir, args, depth_bed_file,gene_mean_depth_dict)

    generate_output_file(all_result_dict, gene_mean_depth_dict, gene_phase_dict, new_result, args["k"], focus_gene_list, variant_num_dict)
    # load_raw_result(raw_result, gene_mean_depth_dict, gene_phase_dict, new_result, args["k"], focus_gene_list)
    print (f"result is in {new_result}")