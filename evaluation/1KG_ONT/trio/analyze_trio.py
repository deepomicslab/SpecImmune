import pandas as pd
import numpy as np
import allel
import scipy
from scipy.spatial import distance

import sys, os
# sys.path.insert(0, sys.path[0]+'/../')
sys.path.insert(0, sys.path[0]+'/../../')
from four_field_compare import convert_field, has_intersection

six_trios = [['2418', 'NA19828', 'NA19818', 'NA19819'],
             ['CLM16', 'HG01258', 'HG01256', 'HG01257'],
             ['1463-Paternal', 'NA12877', 'NA12889', 'NA12890'],
             ['1463-Maternal', 'NA12878', 'NA12891', 'NA12892'],
             ['SH006', 'HG00420', 'HG00418', 'HG00419'],
             ['Y077', 'NA19129', 'NA19128', 'NA19127'],
            ]  # [trio ID, child, father, mother]
            

class My_allele:

    def __init__(self, gene, allele, score, depth, phase_set,variant_num=0, hete_variant_num=0):
        self.gene = gene
        self.allele = allele
        self.score = score
        self.depth = depth
        self.phase_set = phase_set
        self.variant_num = variant_num
        self.hete_variant_num = hete_variant_num


def group_genes(locus):
    if len(locus.split("-")) == 2:
        gene = locus.split("-")[1]
    else:
        gene = locus
    if gene in class1_genes:
        return "Class I"
    if gene in class1_pseudogenes:
        return "Pseudogenes"
    if gene in class2_genes:
        return "Class II"
    if gene in non_hla_genes:
        return "Non-HLA"
    return "other"

## read the file using pandas, sep is \t
def get_super_pop(super_pop_file):
    df = pd.read_csv(super_pop_file, sep='\t')
    super_pop_dict = {}
    for index, row in df.iterrows():
        if  pd.isna(row['Population Code']) or pd.isna(row['Super Population']):
            continue
        super_pop_dict[row['Population Code']] = row['Super Population']
    # print (super_pop_dict)
    return super_pop_dict

def get_sample_pop(sample_pop_file):
    # read excel file into df
    df = pd.read_excel(sample_pop_file)
    sample_pop_dict = {}
    for index, row in df.iterrows():
        sample_pop_dict[row['Sample']] = row['Population']
    # print (sample_pop_dict)
    return sample_pop_dict

def read_alleles(allele_file, super_pop_dict, sample_pop_dict,read_num_cutoff=10,field=8):
    df = pd.read_csv(allele_file)
    alleles_dict = {}
    pop_alleles_dict = {}
    pop_sample_num = {}
    alleles_gene_dict = {}
    alleles_sample_dict = {}

    for index, row in df.iterrows():
        if row['Sample'] not in sample_pop_dict:
            # print ("no pop info for sample", row['Sample'])
            # continue
            pop = 'Non-pop'
            super_pop = 'Non-pop'
        else:
            pop = sample_pop_dict[row['Sample']]
            super_pop = super_pop_dict[pop]
            # print (row['Sample'], sample_pop_dict[row['Sample']])

        if row['Sample'] not in alleles_sample_dict:
            alleles_sample_dict[row['Sample']] = []
        
        if row['Locus'] == "HFE":
            continue
        
        read_num = int(row['Reads_num'])
        if read_num < read_num_cutoff:
            continue
        if super_pop not in alleles_dict:
            alleles_dict[super_pop] = {}
        if row['Locus'] not in alleles_dict[super_pop]:
            alleles_dict[super_pop][row['Locus']] = []
        ## check if empty
        if row['Genotype'] != "NA" and row['Genotype'] != "" and not pd.isna(row['Genotype']):
            genotype = row['Genotype']
            # genotype = genotype.split(";")[0]  ## whether retain all the information
            genotype = convert_field([genotype], field)[0]
            alleles_dict[super_pop][row['Locus']].append(genotype)
            if row['Locus'] not in alleles_gene_dict:
                alleles_gene_dict[row['Locus']] = []
            if pop not in pop_alleles_dict:
                pop_alleles_dict[pop] = {}
            if pop not in pop_sample_num:
                pop_sample_num[pop] = set()
            pop_sample_num[pop].add(row['Sample'])
            if row['Locus'] not in pop_alleles_dict[pop]:
                pop_alleles_dict[pop][row['Locus']] = []
            pop_alleles_dict[pop][row['Locus']].append(genotype)
            alleles_gene_dict[row['Locus']].append(genotype)
            alleles_sample_dict[row['Sample']].append(genotype)
    print ("considered_sample_num", len(alleles_sample_dict), "gene num", len(alleles_gene_dict))
    return alleles_dict, alleles_gene_dict, alleles_sample_dict,pop_alleles_dict,pop_sample_num

def read_alleles_cyp(allele_file, super_pop_dict, sample_pop_dict,read_num_cutoff=10,field=8):
    df = pd.read_csv(allele_file)
    alleles_dict = {}
    pop_alleles_dict = {}
    alleles_gene_dict = {}
    alleles_sample_dict = {}
    sample_phenotype_dict = {}

    for index, row in df.iterrows():
        if row['Sample'] not in sample_pop_dict:
            # print ("no pop info for sample", row['Sample'])
            # continue
            pop = 'Non-pop'
            super_pop = 'Non-pop'
        else:
            pop = sample_pop_dict[row['Sample']]
            super_pop = super_pop_dict[pop]
            # print (row['Sample'], sample_pop_dict[row['Sample']])
        
        # if row['Sample'] == "HG01341":
        #     print ("HG01341 has unkown allele")
        #     continue

        if row['Sample'] not in alleles_sample_dict:
            alleles_sample_dict[row['Sample']] = []
        
        if row['Locus'] == "HFE":
            continue
        
        read_num = int(row['Reads_num'])
        if read_num < read_num_cutoff:
            continue
        if super_pop not in alleles_dict:
            alleles_dict[super_pop] = {}
        if row['Locus'] not in alleles_dict[super_pop]:
            alleles_dict[super_pop][row['Locus']] = []
        ## check if empty
        if row['Genotype'] != "NA" and row['Genotype'] != "" and not pd.isna(row['Genotype']) and row['Genotype'] != "unknown":
            genotype = row['Genotype']
            genotype = genotype.split(";")[0]
            genotype = "CYP2D6"+genotype
            # genotype = genotype.replace("*", "CYP2D6*")
            # genotype = convert_field([genotype], field)[0]
            alleles_dict[super_pop][row['Locus']].append(genotype)
            if row['Locus'] not in alleles_gene_dict:
                alleles_gene_dict[row['Locus']] = []
            if pop not in pop_alleles_dict:
                pop_alleles_dict[pop] = {}
            if row['Locus'] not in pop_alleles_dict[pop]:
                pop_alleles_dict[pop][row['Locus']] = []
            pop_alleles_dict[pop][row['Locus']].append(genotype)
            alleles_gene_dict[row['Locus']].append(genotype)
            alleles_sample_dict[row['Sample']].append(genotype)

            sample_phenotype_dict[row['Sample']] = [row['Phenotype'], row['Activity_score']]

    print ("considered_sample_num", len(sample_phenotype_dict), "gene num", len(alleles_gene_dict))
    # print ("HG00365" in sample_phenotype_dict)
    return alleles_dict, alleles_gene_dict, alleles_sample_dict,pop_alleles_dict,sample_phenotype_dict

def read_alleles_vdj(allele_file, super_pop_dict, sample_pop_dict,read_num_cutoff=10,identity_cutoff=99):
    df = pd.read_csv(allele_file)
    alleles_dict = {}
    pop_alleles_dict = {}
    alleles_gene_dict = {}
    alleles_sample_dict = {}
    sample_gene_dict = {}
    allele_num = 0

    for index, row in df.iterrows():
        if row['Sample'] not in sample_pop_dict:
            # print ("no pop info for sample", row['Sample'])
            # continue
            pop = 'Non-pop'
            super_pop = 'Non-pop'
        else:
            pop = sample_pop_dict[row['Sample']]
            super_pop = super_pop_dict[pop]
        
        read_num = int(row['depth'])
        if read_num < read_num_cutoff:
            continue

        if float(row['score_1']) < identity_cutoff or float(row['score_2']) < identity_cutoff:
            continue

        if super_pop not in alleles_dict:
            alleles_dict[super_pop] = {}
        if row['gene'] not in alleles_dict[super_pop]:
            alleles_dict[super_pop][row['gene']] = []
        ## check if empty
        if row['allele_1'] != "NA" and row['allele_1'] != "" and not pd.isna(row['allele_1']) and row['allele_1'] != "unknown" and\
            row['allele_2'] != "NA" and row['allele_2'] != "" and not pd.isna(row['allele_2']) and row['allele_2'] != "unknown":

            # allele_1 = My_allele(row['gene'], row['allele_1'], row['score_1'], row['depth'], row['phase_set'],row['variant_num'],row['hete_variant_num'])
            # allele_2 = My_allele(row['gene'], row['allele_2'], row['score_2'], row['depth'], row['phase_set'],row['variant_num'],row['hete_variant_num'])
            allele_1 = row['allele_1']
            allele_2 = row['allele_2']
            genotype = [allele_1, allele_2]

            if row['Sample'] not in alleles_sample_dict:
                alleles_sample_dict[row['Sample']] = []
            if row['Sample'] not in sample_gene_dict:
                sample_gene_dict[row['Sample']] = {}

            sample_gene_dict[row['Sample']][row['gene']] = genotype

            alleles_dict[super_pop][row['gene']] = genotype
            if row['gene'] not in alleles_gene_dict:
                alleles_gene_dict[row['gene']] = []
            if pop not in pop_alleles_dict:
                pop_alleles_dict[pop] = {}
            if row['gene'] not in pop_alleles_dict[pop]:
                pop_alleles_dict[pop][row['gene']] = []

            pop_alleles_dict[pop][row['gene']] += genotype
            alleles_gene_dict[row['gene']] += genotype
            alleles_sample_dict[row['Sample']] += genotype

            allele_num += 2

    print ("considered_sample_num", len(alleles_sample_dict), "gene num", len(alleles_gene_dict),"allele num", allele_num)
    return alleles_dict, alleles_gene_dict, alleles_sample_dict,pop_alleles_dict,sample_gene_dict

def split_alleles(alleles):
    gene_result = {}
    for genotypes in alleles:
        geno_list = genotypes.split(";")
        gene = geno_list[0].split("*")[0]
        if gene not in gene_result:
            gene_result[gene] = []
        gene_result[gene].append(geno_list)
    return gene_result

def check_consistent(child_gene_result, father_gene_result, mother_gene_result):
    all_allele_num = 0
    all_consistent_num = 0
    for locus in child_gene_result:
        if len(child_gene_result[locus]) != 2:
            continue
        if locus not in father_gene_result or locus not in mother_gene_result:
            continue

        fir, sec = 0, 0
        if has_intersection(child_gene_result[locus][0], father_gene_result[locus][0] + father_gene_result[locus][1]):
            fir += 1
        if has_intersection(child_gene_result[locus][1], mother_gene_result[locus][0] + mother_gene_result[locus][1]):
            fir += 1   

        
        if has_intersection(child_gene_result[locus][1], father_gene_result[locus][0] + father_gene_result[locus][1]):
            sec += 1
        if has_intersection(child_gene_result[locus][0], mother_gene_result[locus][0] + mother_gene_result[locus][1]):
            sec += 1  
        
        consistent_num = max(fir, sec)
        all_consistent_num += consistent_num
        all_allele_num += 2
    print ("consistent_num", all_consistent_num, "allele_num", all_allele_num)
    return all_consistent_num, all_allele_num

def each_trio(alleles_sample_dict):
    total_all_consistent_num, total_all_allele_num = 0, 0
    for trio in six_trios:
        trio_id = trio[0]
        child = trio[1]
        father = trio[2]
        mother = trio[3]
        # print (alleles_sample_dict[child])
        child_gene_result = split_alleles(alleles_sample_dict[child])
        father_gene_result = split_alleles(alleles_sample_dict[father])
        mother_gene_result = split_alleles(alleles_sample_dict[mother])
        all_consistent_num, all_allele_num = check_consistent(child_gene_result, father_gene_result, mother_gene_result)
        total_all_consistent_num += all_consistent_num
        total_all_allele_num += all_allele_num
        # break
    consistent_freq = total_all_consistent_num/total_all_allele_num
    print ("consistent_freq", consistent_freq, "consistent_num", total_all_consistent_num, "allele_num", total_all_allele_num)

if __name__ == "__main__":
    super_pop_file = "../hla/20131219.populations.tsv"
    sample_pop_file = "../hla/20130606_sample_info.xlsx"

    super_pop_dict = get_super_pop(super_pop_file)
    sample_pop_dict = get_sample_pop(sample_pop_file)

    alleles_dict, alleles_gene_dict_kir, alleles_sample_dict_kir,pop_alleles_dict,pop_sample_num = read_alleles("../kir/merged_samples.csv", super_pop_dict, sample_pop_dict, 10, 8)
    alleles_dict, alleles_gene_dict_hla, alleles_sample_dict_hla,pop_alleles_dict,pop_sample_num = read_alleles("../hla/speclong_res_merged_samples.csv", super_pop_dict, sample_pop_dict, 10, 8)
    alleles_dict, alleles_gene_dict_cyp, alleles_sample_dict_cyp,pop_alleles_dict,sample_phenotype_dict = read_alleles_cyp("../cyp/cyp_1k_all.csv", super_pop_dict, sample_pop_dict, 10, 8)
    alleles_dict, alleles_gene_dict_vdj, alleles_sample_dict_vdj,pop_alleles_dict,sample_gene_dict = read_alleles_vdj("../ig_tr/merged_samples.ig_tr.csv", super_pop_dict, sample_pop_dict, 10, 99)
    # print (alleles_sample_dict_cyp)
    # merge the two gene dict
    alleles_gene_dict_cyp = {"CYP2D6":1}
    alleles_gene_dict = {**alleles_gene_dict_kir, **alleles_gene_dict_hla, **alleles_gene_dict_cyp, **alleles_gene_dict_vdj}

    # alleles_sample_dict = {}
    # for sample in alleles_sample_dict_kir:
    #     if sample in alleles_sample_dict_hla and sample in alleles_sample_dict_kir and sample in alleles_sample_dict_cyp and sample in alleles_sample_dict_vdj:
    #         alleles_sample_dict[sample] = alleles_sample_dict_kir[sample] + alleles_sample_dict_hla[sample]+ alleles_sample_dict_cyp[sample]+ alleles_sample_dict_vdj[sample]

    alleles_sample_dict = {}
    for sample in alleles_sample_dict_hla:
        alleles_sample_dict[sample] = []
        if sample in alleles_sample_dict_hla:
            alleles_sample_dict[sample] += alleles_sample_dict_hla[sample]
        if sample in alleles_sample_dict_kir:
            alleles_sample_dict[sample] += alleles_sample_dict_kir[sample]
        if sample in alleles_sample_dict_cyp:
            alleles_sample_dict[sample] += alleles_sample_dict_cyp[sample]
        if sample in alleles_sample_dict_vdj:
            alleles_sample_dict[sample] += alleles_sample_dict_vdj[sample]
    print ("final sample num", len(alleles_sample_dict))
    each_trio(alleles_sample_dict)