import sys

from four_field_compare import main_TCR,cal_gene_len
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
    tcr_result_dir = "/mnt/d/HLAPro_backup/Nanopore_optimize/vdj_results_tcr/"
    cutoff = 5
    main_TCR(tcr_result_dir, benchmark_result_dir,gene_list,gene_class,cutoff)



