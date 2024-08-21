from collections import defaultdict
import sys, os

sys.path.insert(0, sys.path[0]+'/../../scripts/')

from bed_objects import Bed_db

def load_vdj_result(raw_result, sample, cutoff=0):
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
    phase_set_dict = defaultdict(list)
    # sample = 'test'
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
            phase_set = line[2]
            if allele1 != allele2:
                phase_set_dict[phase_set].append([allele1, allele2])

            if sample not in store_alleles_dict:
                store_alleles_dict[sample] = {}
            if gene not in store_alleles_dict[sample]:
                store_alleles_dict[sample][gene] = []
            store_alleles_dict[sample][gene] = [[allele1], [allele2]]
    
    # for phase_set, alleles in phase_set_dict.items():
    #     print (phase_set, len(alleles))
    return store_alleles_dict

# raw_result = '/mnt/d/HLAPro_backup/Nanopore_optimize/vdj_results2/NA19240/NA19240.IG_TR_typing_result.txt'
# load_vdj_result(raw_result)
bed_db = Bed_db()
hg38_gene_info = bed_db.get_hg38_gene_interval() 
gene_list_sorted_by_pos = list(hg38_gene_info.keys())
focus_gene_list = gene_list_sorted_by_pos[345:367]

for sample in ['NA19238','NA19239','NA19240']:
    raw_result = f'/mnt/d/HLAPro_backup/Nanopore_optimize/vdj_results2/{sample}/{sample}.IG_TR_typing_result.txt'
    store_alleles_dict = load_vdj_result(raw_result,sample)
    for gene in focus_gene_list:
        if gene in store_alleles_dict[sample]:
            print (sample, gene, store_alleles_dict[sample][gene])
        else:
            print (sample, gene, 'NA')
