import sys
import pandas as pd

from four_field_compare import main_TCR,cal_gene_len,main_vdj_hgscv2,main_pacbio,main_kir
sys.path.insert(0, sys.path[0]+'/../scripts/')

from get_lite_db import convert_field_for_allele
from determine_gene import get_focus_gene
from db_objects import My_db

def cal_total_accuracy(sum_result_file):
    ## read the file using pandas
    df = pd.read_csv(sum_result_file)
    # enumerate the df
    depth_dict = {} 
    for index, row in df.iterrows():
        
        if row['cutoff'] not in depth_dict:
            depth_dict[row['cutoff']] = [0, 0, 0]
        depth_dict[row['cutoff']][0] += row['correct']
        depth_dict[row['cutoff']][1] += row['total']
    for cutoff in depth_dict:
        print (cutoff, depth_dict[cutoff], depth_dict[cutoff][0]/depth_dict[cutoff][1])
        depth_dict[cutoff][2] = depth_dict[cutoff][0]/depth_dict[cutoff][1]
    return depth_dict
        

if __name__ == "__main__":

    gene_class = 'KIR'
    step = 2   ### 1 or 2, assess result in step 1 or step 2
    field = 3
    db_dir = f"../db/{gene_class}/"
    gene_list, interval_dict =  get_focus_gene(gene_class)
    gene_mean_len, allele_length_dict = cal_gene_len(db_dir)
    # main_pacbio(gene_list, truth_dir, result_dir, gene_class, step)

    benchmark_result_dir = "kir_results/"
    data = []
    for field in [1,2,3]:
        # truth_dir = "/mnt/d/HLAPro_backup/Nanopore_optimize/hgscv2_truth_bwa/"
        # result_dir = "/mnt/d/HLAPro_backup/Nanopore_optimize/KIR_HGSCV2_hifi2/"
        # sum_result_file = benchmark_result_dir + "HGSCV2_hifi.csv"
        # main_kir(gene_list, truth_dir, result_dir, allele_length_dict, sum_result_file, field)

        # truth_dir = "/mnt/d/HLAPro_backup/Nanopore_optimize/hprc_truth_bwa/"
        # result_dir = "/mnt/d/HLAPro_backup/Nanopore_optimize/KIR_HPRC_hifi/"
        # sum_result_file = benchmark_result_dir + "HPRC_hifi_2.csv"
        # main_kir(gene_list, truth_dir, result_dir, allele_length_dict, sum_result_file, field)

        truth_dir = "/home/wangshuai/00.hla/long/experiments/upload_truth/hprc_truth_bwa_zip/"
        result_dir = "/home/wangshuai/00.hla/long/experiments/kir/kir_results/hprc_hifi2/"
        sum_result_file = benchmark_result_dir + "HPRC_hifi.csv"
        main_kir(gene_list, truth_dir, result_dir, allele_length_dict, sum_result_file, field)

        truth_dir = "/home/wangshuai/00.hla/long/experiments/upload_truth/hprc_truth_bwa_zip/"
        result_dir = "/home/wangshuai/00.hla/long/experiments/kir/kir_results/hprc_ont/"
        sum_result_file = benchmark_result_dir + "HPRC_ont.csv"
        main_kir(gene_list, truth_dir, result_dir, allele_length_dict, sum_result_file, field)

        truth_dir = "/home/wangshuai/00.hla/long/experiments/upload_truth/hgscv2_truth_bwa_zip/"
        result_dir = "/home/wangshuai/00.hla/long/experiments/kir/kir_results/hgscv2_clr/"
        sum_result_file = benchmark_result_dir + "hgscv2_clr.csv"
        main_kir(gene_list, truth_dir, result_dir, allele_length_dict, sum_result_file, field)
 


        # """ 
    #     sum_result_file = benchmark_result_dir + "HPRC_hifi.csv"
    #     result_file = f'{sum_result_file[:-4]}_field{field}_result.csv'
    #     depth_dict = cal_total_accuracy(result_file)
    #     for cutoff in depth_dict:
    #         data.append([cutoff, depth_dict[cutoff][0], depth_dict[cutoff][1], depth_dict[cutoff][2], 'HPRC HiFi', field])


    #     sum_result_file = benchmark_result_dir + "HGSCV2_hifi.csv"
    #     result_file = f'{sum_result_file[:-4]}_field{field}_result.csv'
    #     depth_dict = cal_total_accuracy(result_file)
    #     for cutoff in depth_dict:
    #         data.append([cutoff, depth_dict[cutoff][0], depth_dict[cutoff][1], depth_dict[cutoff][2], 'HGSCV2 HiFi', field])

    #     sum_result_file = benchmark_result_dir + "HPRC_ont.csv"
    #     result_file = f'{sum_result_file[:-4]}_field{field}_result.csv'
    #     depth_dict = cal_total_accuracy(result_file)
    #     for cutoff in depth_dict:
    #         data.append([cutoff, depth_dict[cutoff][0], depth_dict[cutoff][1], depth_dict[cutoff][2], 'HPRC ONT', field])

    #     sum_result_file = benchmark_result_dir + "hgscv2_clr.csv"
    #     result_file = f'{sum_result_file[:-4]}_field{field}_result.csv'
    #     depth_dict = cal_total_accuracy(result_file)
    #     for cutoff in depth_dict:
    #         data.append([cutoff, depth_dict[cutoff][0], depth_dict[cutoff][1], depth_dict[cutoff][2], 'HGSCV2 CLR', field])
    #     # """

    # # """
    # df = pd.DataFrame(data, columns = ['depth', 'correct', 'total', 'accuracy', 'dataset', 'field'])
    # df.to_csv(benchmark_result_dir + "/all_loci_depth.csv", index=False)
    # # """
