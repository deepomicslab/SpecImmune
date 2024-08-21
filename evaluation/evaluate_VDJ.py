import sys

from four_field_compare import main_TCR,cal_gene_len,main_vdj_hgscv2
sys.path.insert(0, sys.path[0]+'/../scripts/')

from get_lite_db import convert_field_for_allele
from determine_gene import get_focus_gene
from db_objects import My_db


if __name__ == "__main__":

    gene_class = 'IG_TR'
    step = 2   ### 1 or 2, assess result in step 1 or step 2
    db_dir = f"../db/{gene_class}/"
    gene_list, interval_dict =  get_focus_gene(gene_class)
    gene_mean_len, allele_length_dict = cal_gene_len(db_dir)
    # main_pacbio(gene_list, truth_dir, result_dir, gene_class, step)

    benchmark_result_dir = "vdj_results/"

    # sum_result_file = benchmark_result_dir + "tcr_11.csv"
    # tcr_result_dir = "/mnt/d/HLAPro_backup/Nanopore_optimize/vdj_results_tcr4/"
    # cutoff = 2
    # main_TCR(tcr_result_dir, benchmark_result_dir,gene_list,gene_class,sum_result_file,cutoff)

    #### evaluation in HGSCV2
    gene_class = "IG_TR"
    truth_dir = "/mnt/d/HLAPro_backup/Nanopore_optimize/hgscv2_truth_bwa/"
    result_dir = "/mnt/d/HLAPro_backup/Nanopore_optimize/vdj_results2/"
    sum_result_file = benchmark_result_dir + "HGSCV2_hifi.csv"
    # cutoff = 2
    main_vdj_hgscv2(gene_list, truth_dir, result_dir, allele_length_dict, sum_result_file, gene_class)



