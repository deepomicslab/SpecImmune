
import sys
import argparse
import pandas as pd

sys.path.insert(0, sys.path[0]+'/../scripts/')

from four_field_compare import assess_sim_module_hprc_pangenome
from determine_gene import get_focus_gene
from trans_skirt2csv import parse_kir_csv, save_results



args = {}
args['i'] = 'KIR'
gene_list, interval_dict =  get_focus_gene(args)
truth_dir="/data6/wangxuedong/add_experiment_HPRC_KIR/truth_alleles_per_sample/"
res_dir="/data6/wangxuedong/add_experiment_HPRC_KIR/out_dir/"
# samples=HG002 HG00438 HG005 HG00621 HG00673 HG00733 HG00735 HG00741 HG01071 HG01106 HG01123 HG01175 HG01258 HG01358 HG01361 HG01891 HG01928 HG01952 HG01978 HG02148 HG02257 HG02486 HG02559 HG02572 HG02622 HG02630 HG02717 HG02886 HG03453 HG03471 HG03516 HG03540 HG03579 HG01109 HG01243 HG02055 HG02080 HG02109 HG02145 HG02723 HG02818 HG03098 HG03486 HG03492 NA18906 NA19240 NA20129 NA21309"
samples=['HG00438', 'HG00621', 'HG00673', 'HG00733', 'HG00735', 'HG00741', 'HG01071', 'HG01106', 'HG01123', 'HG01175', 'HG01258', 'HG01358', 'HG01361', 'HG01891', 'HG01928', 'HG01952', 'HG01978', 'HG02148', 'HG02257', 'HG02486', 'HG02559', 'HG02572', 'HG02622', 'HG02630', 'HG02717', 'HG02886', 'HG03453', 'HG03516', 'HG03540', 'HG03579', 'HG01109', 'HG01243', 'HG02055', 'HG02080', 'HG02109', 'HG02145', 'HG02723', 'HG02818', 'HG03098', 'HG03486', 'HG03492', 'NA18906', 'NA19240', 'NA20129', 'NA21309']
print("we have ", len(samples), " samples")

depth_threshold = 5
depth_file="/data6/wangxuedong/add_experiment_HPRC_KIR/genes_depth.csv"
# sample,gene,haplotype,average_depth
# HG00438,KIR2DL1,0,29.09416
# HG00438,KIR2DL1,1,8.49839
# HG00438,KIR2DL3,0,22.19242
# HG00438,KIR2DL3,1,10.01782
# HG00438,KIR2DL4,0,25.87555
# HG00438,KIR2DL4,1,14.51364
depth_df = pd.read_csv(depth_file)

focus_genes=["KIR3DL1","KIR3DL2","KIR3DS1"]



gene_count_table = {}
for sample in samples:

    true=f'{truth_dir}/{sample}.csv'
    infer= f'{res_dir}/{sample}/{sample}/out/merged.hap.csv'
    infer_f = f'{res_dir}/{sample}/{sample}/out/merged.hap.format.csv'
    infer_f2 = f'{res_dir}/{sample}/{sample}/{sample}.KIR.final.type.result.formatted.txt'
    results = parse_kir_csv(infer)
    save_results(results, infer_f)

    spec_gene_accuracy_dict = assess_sim_module_hprc_pangenome(true, infer_f, gene_list, args['i'], '',1)
    # spec_gene_accuracy_dict = assess_sim_module_hprc_pangenome(true, infer_f2, gene_list, args['i'], '',3)
    for gene in spec_gene_accuracy_dict:
        if gene not in focus_genes:
            continue
        print (sample, gene, spec_gene_accuracy_dict[gene])
        # get the depth for haplotype 0 of gene of sample, if it exists, else 0
        depth = depth_df[(depth_df['sample'] == sample) & (depth_df['gene'] == gene) & (depth_df['haplotype'] == '0')]
        if len(depth) > 0:
            depth = float(depth.iloc[0]['average_depth'])
        else:
            depth = 0
        # get the depth for haplotype 1 of gene of sample, if it exists, else 0
        depth1 = depth_df[(depth_df['sample'] == sample) & (depth_df['gene'] == gene) & (depth_df['haplotype'] == '1')]
        if len(depth1) > 0:
            depth1 = float(depth1.iloc[0]['average_depth'])
        else:
            depth1 = 0
        all_depth= depth + depth1
        if all_depth < depth_threshold:
            continue


        if gene not in gene_count_table:
            gene_count_table[gene] = [0,0]
        gene_count_table[gene][0] += spec_gene_accuracy_dict[gene][1]
        gene_count_table[gene][1] += spec_gene_accuracy_dict[gene][2]

print(gene_count_table)
all_match=0
all_count=0
# calculate the average
for gene in gene_count_table:
    acc = gene_count_table[gene][0] /  gene_count_table[gene][1]
    print (gene, acc)
    all_match += gene_count_table[gene][0]
    all_count += gene_count_table[gene][1]
print ("all match", all_match)
print ("all count", all_count)
print ("all accuracy", all_match / all_count)
