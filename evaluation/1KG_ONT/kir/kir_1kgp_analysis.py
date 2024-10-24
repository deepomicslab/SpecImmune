import pandas as pd
import numpy as np
import allel
import scipy
from scipy.spatial import distance

import sys, os
# sys.path.insert(0, sys.path[0]+'/../')
sys.path.insert(0, sys.path[0]+'/../../')
from four_field_compare import convert_field

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
            genotype = genotype.split(";")[0]
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


def count_freq(list, count='no'):
    ## count the frequency of each allele
    freq_dict = {}
    total = 0
    for allele in list:
        if allele not in freq_dict:
            freq_dict[allele] = 0
        total += 1
        freq_dict[allele] += 1
    if count == 'yes':
        return freq_dict
    for allele in freq_dict:
        freq_dict[allele] = freq_dict[allele]/total
    return freq_dict

def count_alleles(alleles_dict, freq_file):
    ## count the frequency of each allele
    
    freq_dict = {}
    super_pop_freq_dist = {}
    all_pop_alleles = []
    super_pop_list = list(alleles_dict.keys())
    for super_pop in alleles_dict:
        if super_pop not in super_pop_freq_dist:
            super_pop_freq_dist[super_pop] = {}
        all_alleles = []
        for locus in alleles_dict[super_pop]:
            all_alleles += alleles_dict[super_pop][locus]
            locus_freq_dict = count_freq(alleles_dict[super_pop][locus])
            print (locus_freq_dict)
            super_pop_freq_dist[super_pop].update(locus_freq_dict)
        all_pop_alleles += all_alleles
        # freq_dict = count_freq(all_alleles)
        # super_pop_freq_dist[super_pop] = freq_dict

        # for allele in freq_dict:
        #     data.append([super_pop, allele, freq_dict[allele], str(allele).split("*")[0]])
    uniq_alleles = list(set(all_pop_alleles))
    data = []
    for allele in uniq_alleles:
        pop_freq_list = []
        for super_pop in super_pop_list:
            if allele in super_pop_freq_dist[super_pop]:
                pop_freq_list.append(super_pop_freq_dist[super_pop][allele])
            else:
                pop_freq_list.append(0)
        data.append([allele, str(allele).split("*")[0]] + pop_freq_list)
    df = pd.DataFrame(data, columns=['Allele', 'locus'] + super_pop_list)
    df.to_csv(freq_file, index=False)

def for_histogram(count_alleles, histogram_file):
    data = []
    all_alleles = []
    for super_pop in alleles_dict:
        # if super_pop != 'EAS':
        #     continue
        for locus in alleles_dict[super_pop]:
            all_alleles += alleles_dict[super_pop][locus]
    allele_count = count_freq(all_alleles, count='yes')
    for allele in allele_count:
        data.append([allele, allele_count[allele]])
    df = pd.DataFrame(data, columns=['Allele', 'Count'])
    df.to_csv(histogram_file, index=False)

def MAF_analysis(alleles_gene_dict, MAF_file, gene_allele_file):
    data = []
    gene_count = []
    summary = [0, 0, 0, 0]
    allele_alleles = []
    for gene in alleles_gene_dict:
        alleles = alleles_gene_dict[gene]
        gene_count.append([gene, len(set(alleles)), 'uniq'])
        gene_count.append([gene, len(alleles), 'all'])
        allele_alleles += alleles
        ## get the unique alleles, and count the frequency, count the allele number with frequency < 0.05, and the allele number with fequency >= 0.05 
        allele_count = count_freq(alleles, count='no')
        allele_count_low = 0
        allele_count_high = 0
        allele_count_rare = 0
        this_allele_num = 0
        for allele in allele_count:
            if allele_count[allele] > 0.05:
                allele_count_high += 1
                summary[0] += 1
            elif allele_count[allele] > 0.005:
                allele_count_low += 1
                summary[1] += 1
            else:
                allele_count_rare += 1
                summary[2] += 1
            this_allele_num += 1
            summary[3] += 1
        gene_group = gene[:4]
        # data.append([gene, len(alleles), len(allele_count), allele_count_low, allele_count_high])
        data.append([gene, gene_group, allele_count_rare, allele_count_rare/this_allele_num, 'rare'])
        data.append([gene, gene_group, allele_count_low, allele_count_low/this_allele_num,'low-frequency'])
        data.append([gene, gene_group, allele_count_high, allele_count_high/this_allele_num,'common'])

    print (len(set(allele_alleles)))
    df = pd.DataFrame(data, columns=['Gene', 'Gene Group', 'Count', 'Frequency', 'Group'])
    df.to_csv(MAF_file, index=False)
    print (summary, 'high low rare', np.array(summary)/summary[3])

    df = pd.DataFrame(gene_count, columns=['Gene', 'Count', 'Type'])
    df.to_csv(gene_allele_file, index=False)

def cumulative_analysis(alleles_sample_dict, super_pop_dict, cumulative_file):
    all_uniq_alleles = []
    num = 0
    data = []
    for sample in alleles_sample_dict:
        all_uniq_alleles += alleles_sample_dict[sample]
        num += 1
        data.append([num, len(set(all_uniq_alleles))])
    df = pd.DataFrame(data, columns=['Sample_Num', 'Allele_Num'])
    df.to_csv(cumulative_file, index=False)
    print (len(set(all_uniq_alleles)))

def calculate_fst(allele_freq_group_a, allele_freq_group_b):
    p_bar = (np.array(allele_freq_group_a) + np.array(allele_freq_group_b)) / 2
    p_i_group_a = np.mean(allele_freq_group_a)
    p_i_group_b = np.mean(allele_freq_group_b)
    
    Hs = np.mean((np.array(allele_freq_group_a) - p_i_group_a)**2 + (np.array(allele_freq_group_b) - p_i_group_b)**2)
    Ht = np.mean((p_bar - p_i_group_a)**2 + (p_bar - p_i_group_b)**2)
    
    fst = (Ht - Hs) / Ht
    
    return fst

def Fst_analysis(alleles_dict, fst_file):
    super_pop_freq_dist = {}
    all_pop_alleles = []
    super_pop_list = list(alleles_dict.keys())
    for super_pop in alleles_dict:
        if super_pop == 'Non-pop':
            continue
        if super_pop not in super_pop_freq_dist:
            super_pop_freq_dist[super_pop] = {}
        all_alleles = []
        for locus in alleles_dict[super_pop]:
            all_alleles += alleles_dict[super_pop][locus]
            all_pop_alleles += alleles_dict[super_pop][locus]

        locus_freq_dict = count_freq(all_alleles)
        super_pop_freq_dist[super_pop].update(locus_freq_dict)
    
    uniq_alleles = list(set(all_pop_alleles))
    allele_freq  = []
    group_names = []
    for super_pop in super_pop_freq_dist:
        freq = []
        for allele in uniq_alleles:
            if allele not in super_pop_freq_dist[super_pop]:
                freq.append(0) 
            else:
                freq.append(super_pop_freq_dist[super_pop][allele])
        allele_freq.append(freq)
        group_names.append(super_pop)

    # Calculate pairwise Fst values using scikit-allel
    # pairwise_fst_values = calculate_pairwise_fst(allele_freq)

    data = []
    # Print the pairwise Fst values for each pair of groups
    for i in range(len(group_names)):
        for j in range(len(group_names)):
            fst_value = calculate_fst(allele_freq[i], allele_freq[j])
            ## calculate correlation for the two groups using scipy
            result = scipy.stats.pearsonr(allele_freq[i], allele_freq[j])
            jsd = distance.jensenshannon(allele_freq[i], allele_freq[j])
            print(f"Fst value between {group_names[i]} and {group_names[j]}: {jsd}")
            if super_pop_dict[group_names[i]] == super_pop_dict[group_names[j]]:
                compare = 'intra'
            else:
                compare = 'inter'
            data.append([group_names[i], group_names[j], jsd, compare])
    df = pd.DataFrame(data, columns=['Group1', 'Group2', 'JSD', 'compare'])
    df.to_csv(fst_file, index=False)

def sort_pop(super_pop_dict, color_file):
    pop_dict = {}
    for pop in super_pop_dict:
        if super_pop_dict[pop] not in pop_dict:
            pop_dict[super_pop_dict[pop]] = []
        pop_dict[super_pop_dict[pop]].append(pop)
    pop_list = []
    for super_pop in pop_dict:
        print (super_pop, pop_dict[super_pop])
        pop_dict[super_pop] = sorted(pop_dict[super_pop])
        pop_list += pop_dict[super_pop]
    print (pop_list)
    data = []
    for pop in pop_list:
        data.append([pop, 1, super_pop_dict[pop]])
    df = pd.DataFrame(data, columns=['y', 'x', 'group'])
    df.to_csv(color_file, index=False)
    return pop_list
            
def pop_occur_freq(pop_alleles_dict,pop_sample_num):
    data = []
    for super_pop in pop_alleles_dict:
        if super_pop == 'Non-pop':
            continue

        for locus in pop_alleles_dict[super_pop]:
            count =  len(pop_alleles_dict[super_pop][locus])
            freq= count/(len(pop_sample_num[super_pop])*2)
            data.append([super_pop, locus, count, freq])
            # print (super_pop, locus, count, freq)
    df = pd.DataFrame(data, columns=['Group', 'Gene', 'Count','Freq'])
    df.to_csv("pop_occur_freq.csv", index=False)

def LD_analysis(alleles_gene_dict, alleles_sample_dict,LD_dir):
    loci_list =list(alleles_gene_dict.keys())
    for i in range(len(loci_list)):
        for j in range(i+1, len(loci_list)):
            gene1 = loci_list[i]
            gene2 = loci_list[j]

            data = []
            for sample in alleles_sample_dict:
                gene1_alleles = []
                gene2_alleles = []
                for allele in alleles_sample_dict[sample]:
                    if allele.split("*")[0] == gene1:
                        gene1_alleles.append(allele)
                    if allele.split("*")[0] == gene2:
                        gene2_alleles.append(allele)
                if len(gene1_alleles) == 2 and len(gene2_alleles) == 2:
                    # print (gene1_alleles[0], gene1_alleles[1], gene2_alleles[0], gene2_alleles[1])
                    data.append([gene1_alleles[0], gene1_alleles[1], gene2_alleles[0], gene2_alleles[1]])
            if len(data) > 0:
                df = pd.DataFrame(data, columns=[gene1, gene1, gene2, gene2])
                df.to_csv(f"{LD_dir}/{gene1}_{gene2}_LD.csv", index=False)

super_pop_file = "../hla/20131219.populations.tsv"
sample_pop_file = "../hla/20130606_sample_info.xlsx"

allele_file = "merged_samples.csv"
freq_file = "./kir_freq.csv"
histogram_file = "./kir_hist.csv"
gene_allele_file = "./kir_gene_allele.csv"
MAF_file = "./kir_maf.csv"
cumulative_file = "./kir_cumulative.csv"
fst_file = "./kir_fst.csv"
color_file = "./kir_color.csv"
super_pop_dict = get_super_pop(super_pop_file)
sample_pop_dict = get_sample_pop(sample_pop_file)
alleles_dict, alleles_gene_dict, alleles_sample_dict,pop_alleles_dict,pop_sample_num = read_alleles(allele_file, super_pop_dict, sample_pop_dict, 10, 8)


# count_alleles(alleles_dict, freq_file)

# for_histogram(alleles_gene_dict, histogram_file)
# os.system("Rscript plot_hitogram.R")

# MAF_analysis(alleles_gene_dict, MAF_file, gene_allele_file)
# pop_occur_freq(pop_alleles_dict,pop_sample_num)

# cumulative_analysis(alleles_sample_dict, super_pop_dict, cumulative_file)

# Fst_analysis(pop_alleles_dict, fst_file)
# sort_pop(super_pop_dict, color_file)

# LD_dir = "/mnt/d/HLAPro_backup/Nanopore_optimize/1kgp_analysis/kir_LD"
# LD_analysis(alleles_gene_dict, alleles_sample_dict,LD_dir)

# alleles_dict, alleles_gene_dict, alleles_sample_dict,pop_alleles_dict,pop_sample_num = read_alleles("../hla/speclong_res_merged_samples.csv", super_pop_dict, sample_pop_dict, 10, 8)
# LD_dir = "/mnt/d/HLAPro_backup/Nanopore_optimize/1kgp_analysis/hla_LD"
# LD_analysis(alleles_gene_dict, alleles_sample_dict,LD_dir)



###### HLA &&& KIR
"""
alleles_dict, alleles_gene_dict_kir, alleles_sample_dict_kir,pop_alleles_dict,pop_sample_num = read_alleles(allele_file, super_pop_dict, sample_pop_dict, 10, 8)
alleles_dict, alleles_gene_dict_hla, alleles_sample_dict_hla,pop_alleles_dict,pop_sample_num = read_alleles("../hla/speclong_res_merged_samples.csv", super_pop_dict, sample_pop_dict, 10, 8)

# merge the two gene dict
alleles_gene_dict = {**alleles_gene_dict_kir, **alleles_gene_dict_hla}
alleles_sample_dict = {}
for sample in alleles_sample_dict_kir:
    if sample in alleles_sample_dict_hla and sample in alleles_sample_dict_kir:
        alleles_sample_dict[sample] = alleles_sample_dict_kir[sample] + alleles_sample_dict_hla[sample]
LD_dir = "/mnt/d/HLAPro_backup/Nanopore_optimize/1kgp_analysis/hla_kir_LD"
LD_analysis(alleles_gene_dict, alleles_sample_dict,LD_dir)

# data = []
# for gene in alleles_gene_dict_kir:
#     data.append([gene, 1, "KIR"])
# for gene in alleles_gene_dict_hla:
#     data.append([gene, 1, "HLA"])
# df = pd.DataFrame(data, columns=['y', 'x', 'Type'])
# df.to_csv("gene_type.csv", index=False)

"""

# ##### HLA CYP KIR
# alleles_dict, alleles_gene_dict_kir, alleles_sample_dict_kir,pop_alleles_dict,pop_sample_num = read_alleles(allele_file, super_pop_dict, sample_pop_dict, 10, 8)
# alleles_dict, alleles_gene_dict_hla, alleles_sample_dict_hla,pop_alleles_dict,pop_sample_num = read_alleles("../hla/speclong_res_merged_samples.csv", super_pop_dict, sample_pop_dict, 10, 8)
# alleles_dict, alleles_gene_dict_cyp, alleles_sample_dict_cyp,pop_alleles_dict,sample_phenotype_dict = read_alleles_cyp("../cyp/cyp_1k_all.csv", super_pop_dict, sample_pop_dict, 10, 8)

# # print (alleles_sample_dict_cyp)
# # merge the two gene dict
# alleles_gene_dict_cyp = {"CYP2D6":1}
# alleles_gene_dict = {**alleles_gene_dict_kir, **alleles_gene_dict_hla, **alleles_gene_dict_cyp}
# alleles_sample_dict = {}
# for sample in alleles_gene_dict_hla:
#     alleles_sample_dict[sample] = []
#     if sample in alleles_sample_dict_hla:
#         alleles_sample_dict[sample] += alleles_sample_dict_hla[sample]
#     if sample in alleles_sample_dict_kir:
#         alleles_sample_dict[sample] += alleles_sample_dict_kir[sample]
#     if sample in alleles_sample_dict_cyp:
#         alleles_sample_dict[sample] += alleles_sample_dict_cyp[sample]

#     # if sample in alleles_sample_dict_hla and sample in alleles_sample_dict_kir and sample in alleles_sample_dict_cyp:
#     #     alleles_sample_dict[sample] = alleles_sample_dict_kir[sample] + alleles_sample_dict_hla[sample]+ alleles_sample_dict_cyp[sample]
# LD_dir = "/mnt/d/HLAPro_backup/Nanopore_optimize/1kgp_analysis/hla_kir_cyp_LD"
# LD_analysis(alleles_gene_dict, alleles_sample_dict,LD_dir)

# data = [['CYP2D6', 1, "CYP"]]
# for gene in alleles_gene_dict_kir:
#     data.append([gene, 1, "KIR"])
# for gene in alleles_gene_dict_hla:
#     data.append([gene, 1, "HLA"])
# df = pd.DataFrame(data, columns=['y', 'x', 'Type'])
# df.to_csv("cyp_gene_type.csv", index=False)



# # """
# ##### HLA CYP KIR VDJ
alleles_dict, alleles_gene_dict_kir, alleles_sample_dict_kir,pop_alleles_dict,pop_sample_num = read_alleles(allele_file, super_pop_dict, sample_pop_dict, 10, 8)
alleles_dict, alleles_gene_dict_hla, alleles_sample_dict_hla,pop_alleles_dict,pop_sample_num = read_alleles("../hla/speclong_res_merged_samples.csv", super_pop_dict, sample_pop_dict, 10, 8)
alleles_dict, alleles_gene_dict_cyp, alleles_sample_dict_cyp,pop_alleles_dict,sample_phenotype_dict = read_alleles_cyp("../cyp/cyp_1k_all.csv", super_pop_dict, sample_pop_dict, 10, 8)
alleles_dict, alleles_gene_dict_vdj, alleles_sample_dict_vdj,pop_alleles_dict,sample_gene_dict = read_alleles_vdj("../ig_tr/merged_samples.ig_tr.csv", super_pop_dict, sample_pop_dict, 10, 99)
# print (alleles_sample_dict_cyp)
# merge the two gene dict
alleles_gene_dict_cyp = {"CYP2D6":1}
alleles_gene_dict = {**alleles_gene_dict_kir, **alleles_gene_dict_hla, **alleles_gene_dict_cyp, **alleles_gene_dict_vdj}

# # alleles_sample_dict = {}
# # for sample in alleles_sample_dict_kir:
# #     if sample in alleles_sample_dict_hla and sample in alleles_sample_dict_kir and sample in alleles_sample_dict_cyp and sample in alleles_sample_dict_vdj:
# #         alleles_sample_dict[sample] = alleles_sample_dict_kir[sample] + alleles_sample_dict_hla[sample]+ alleles_sample_dict_cyp[sample]+ alleles_sample_dict_vdj[sample]

# alleles_sample_dict = {}
# for sample in alleles_sample_dict_hla:
#     alleles_sample_dict[sample] = []
#     if sample in alleles_sample_dict_hla:
#         alleles_sample_dict[sample] += alleles_sample_dict_hla[sample]
#     if sample in alleles_sample_dict_kir:
#         alleles_sample_dict[sample] += alleles_sample_dict_kir[sample]
#     if sample in alleles_sample_dict_cyp:
#         alleles_sample_dict[sample] += alleles_sample_dict_cyp[sample]
#     if sample in alleles_sample_dict_vdj:
#         alleles_sample_dict[sample] += alleles_sample_dict_vdj[sample]
# print ("final sample num", len(alleles_sample_dict))

####### output the type of all loci of each sample
# sample_type_result = '../chromoMap/sample_all_loci_type.csv'
# f_out = open(sample_type_result, 'w')
# for sample in alleles_sample_dict:
#     sample_alleles = [sample] + alleles_sample_dict[sample]
#     # print (sample, alleles_sample_dict[sample])
#     print (",".join(sample_alleles), file=f_out)
# f_out.close()

# LD_dir = "/mnt/d/HLAPro_backup/Nanopore_optimize/1kgp_analysis/hla_kir_cyp_vdj_LD"
# LD_analysis(alleles_gene_dict, alleles_sample_dict,LD_dir)

data = [['CYP2D6', 1, "CYP", "CYP", "#d9e6eb"]]
for gene in alleles_gene_dict_kir:
    data.append([gene, 1, "KIR", "KIR", "#8f96bd"])
for gene in alleles_gene_dict_hla:
    data.append([gene, 1, "HLA", "HLA", "#9fc3d5"])
for gene in alleles_gene_dict_vdj:
    if gene[:2] == 'IG':
        data.append([gene, 1, "VDJ", 'IG', "#2a347a"])
    else:
        data.append([gene, 1, "VDJ", 'TCR', "#d6d69b"])
df = pd.DataFrame(data, columns=['y', 'x', 'Type', 'subtype', 'color'])
df.to_csv("vdj_gene_type.csv", index=False)
# # """

colors = {
    'HLA': '#9fc3d5',
    'KIR': '#8f96bd',
    'CYP': '#d9e6eb',
    'IG': "#2a347a",
    'TCR': "#d6d69b"}
