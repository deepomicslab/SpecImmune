import pandas as pd
import numpy as np
import allel
import scipy
from scipy.spatial import distance
from collections import defaultdict

import sys, os
from statsmodels.stats.multitest import multipletests
sys.path.insert(0, sys.path[0]+'/../')
# from four_field_compare import convert_field

class1_genes = ['A', 'B', 'C', 'E', 'F', 'G']
class1_pseudogenes = ['H', 'J', 'K', 'L', 'N', 'P', 'S', 'T', 'U', 'V', 'W', 'Y']
class2_genes = ['DRA', 'DQA1', 'DQA2', 'DQB1', 'DQB2', 'DPA1', 'DPA2', 'DPB1', 'DPB2', 'DMA', 'DMB', 'DOA', 'DOB', 
                'DRB1', 'DRB3', 'DRB4', 'DRB5']
non_hla_genes = ['MICA', 'MICB', 'TAP1', 'TAP2']

def get_gene_list():
    HLA = class1_genes + class2_genes + class1_pseudogenes
    new_HLA = ['HLA-' + gene for gene in HLA]
    return new_HLA + non_hla_genes

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

class My_allele:

    def __init__(self, gene, allele, score, depth, phase_set,variant_num=0, hete_variant_num=0):
        self.gene = gene
        self.allele = allele
        self.score = score
        self.depth = depth
        self.phase_set = phase_set
        self.variant_num = variant_num
        self.hete_variant_num = hete_variant_num



def read_alleles(allele_file, super_pop_dict, sample_pop_dict,read_num_cutoff=10,identity_cutoff=99):
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

            allele_1 = My_allele(row['gene'], row['allele_1'], row['score_1'], row['depth'], row['phase_set'],row['variant_num'],row['hete_variant_num'])
            allele_2 = My_allele(row['gene'], row['allele_2'], row['score_2'], row['depth'], row['phase_set'],row['variant_num'],row['hete_variant_num'])
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
    gene_list = get_gene_list()
    for gene in alleles_gene_dict:
        if gene not in gene_list:
            print (gene, "not in the gene_list")
    print (len(gene_list), gene_list)
    summary = [0, 0, 0, 0]
    allele_alleles = []
    for gene in alleles_gene_dict:
    # for gene in gene_list:
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
        gene_group = group_genes(gene)
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
            data.append([group_names[i], group_names[j], jsd])
    df = pd.DataFrame(data, columns=['Group1', 'Group2', 'JSD'])
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
            
def get_allele_proportion(super_pop_dict, sample_pop_dict, alleles_sample_dict,hap_freq_file,hap_count_file):
    data = []
    pop_allele_dict = {}
    hap_count_dict = {}
    all_haps = []
    for sample in alleles_sample_dict:
        if sample not in sample_pop_dict:
            continue
        pop = sample_pop_dict[sample]
        super_pop = super_pop_dict[pop]
        if pop not in pop_allele_dict:
            pop_allele_dict[pop] = []
        pop_allele_dict[pop] += alleles_sample_dict[sample]
        all_haps += alleles_sample_dict[sample]
        for hap in alleles_sample_dict[sample]:
            if hap not in hap_count_dict:
                hap_count_dict[hap] = 0
            hap_count_dict[hap] += 1
    uniq_haps = list(set(all_haps))
    for pop in pop_allele_dict:
        total = len(pop_allele_dict[pop])
        allele_freq = count_freq(pop_allele_dict[pop], "yes")
        for allele in uniq_haps:
            if allele in allele_freq:
                freq = allele_freq[allele]/total
                data.append([pop, allele, allele_freq[allele], total, freq])
            else:
                data.append([pop, allele, 0, total, 0])
    print ("unique haplotypes number:", len(uniq_haps))
    # print (data)
    df = pd.DataFrame(data, columns=['Pop', 'Allele', 'Count', 'Total', 'Frequency'])
    df.to_csv(hap_freq_file, index=False)

    ## sort hap by count using hap_count_dict
    sorted_haps = sorted(hap_count_dict.items(), key=lambda x: x[1], reverse=True)
    print (sorted_haps[:5])

    data = []
    for hap in hap_count_dict:
        data.append([hap, hap_count_dict[hap]])
    df = pd.DataFrame(data, columns=['Hap', 'Count'])
    df.to_csv(hap_count_file, index=False)


def count_vdj(sample_gene_dict, identity_cutoff=99):
    data = []
    
    gene_list = set()
    hete_count_dict = defaultdict(list)
    hete_freq_dict_pop = {}
    hete_group_count_dict = defaultdict(list)
    hete_count_dict_pop = defaultdict(list)
    for sample in sample_gene_dict:

        if sample not in sample_pop_dict:
            pop = 'Non-pop'
            super_pop = 'Non-pop'
        else:
            pop = sample_pop_dict[sample]
            super_pop = super_pop_dict[pop]
        
        if super_pop != 'AFR':
            new_group = 'other'
        else:
            new_group = 'AFR'

        # new_group = pop
        if new_group not in hete_count_dict_pop:
            hete_count_dict_pop[new_group] = {}
        if pop not in hete_freq_dict_pop:
            hete_freq_dict_pop[pop] = {}


        for gene in sample_gene_dict[sample]:
            if sample_gene_dict[sample][gene][0].score < identity_cutoff or sample_gene_dict[sample][gene][1].score < identity_cutoff:
                continue
            gene_list.add(gene)
            group = gene[:3]
            if gene not in hete_count_dict:
                hete_count_dict[gene] = [0, 0]
            if group not in hete_group_count_dict:
                hete_group_count_dict[group] = [0, 0]
            if gene not in hete_count_dict_pop[new_group]:
                hete_count_dict_pop[new_group][gene] = [0, 0]
            if gene not in hete_freq_dict_pop[pop]:
                hete_freq_dict_pop[pop][gene] = [0, 0]

            if sample_gene_dict[sample][gene][0].allele != sample_gene_dict[sample][gene][1].allele:
                hete_count_dict[gene][0] += 1
                hete_group_count_dict[group][0] += 1
                hete_count_dict_pop[new_group][gene][0] += 1
                hete_freq_dict_pop[pop][gene][0] += 1
                allele_num = 2
            else:
                hete_count_dict[gene][1] += 1
                hete_group_count_dict[group][1] += 1
                hete_count_dict_pop[new_group][gene][1] += 1
                hete_freq_dict_pop[pop][gene][1] += 1
                allele_num = 1
            
            data.append([sample, pop, super_pop, gene, group, \
                         sample_gene_dict[sample][gene][0].hete_variant_num,sample_gene_dict[sample][gene][0].variant_num,new_group,allele_num])
        # break


    homo_gene_num = 0
    hete_genes = []
    data3 = []
    for gene in hete_count_dict:
        if hete_count_dict[gene][0] == 0:
            homo_gene_num += 1
        hete_ratio = hete_count_dict[gene][0]/(hete_count_dict[gene][0] + hete_count_dict[gene][1])
        data3.append([gene, gene[:3],gene[:2], hete_count_dict[gene][0], hete_count_dict[gene][1], hete_ratio])
        if hete_ratio > 0.1:
            hete_genes.append(gene)
        # print (gene, hete_count_dict[gene], hete_count_dict[gene][0]/(hete_count_dict[gene][0] + hete_count_dict[gene][1]))
    print (len(hete_count_dict), homo_gene_num, homo_gene_num/len(hete_count_dict), 1-homo_gene_num/len(hete_count_dict))
    df = pd.DataFrame(data3, columns=['Gene', 'Gene_Group', 'Gene_Group2', 'Hete','Homo','Hete_freq'])
    df.to_csv("vdj_hete_freq.csv", index=False)

    # data = [x for x in data if x[3] in hete_genes]
    df = pd.DataFrame(data, columns=['Sample', 'Pop', 'Super Pop', 'Gene', 'Group', 'Hete_Variant_Num','Variant_Num','new_pop','allele_num'])
    df.to_csv("vdj_snp_count.csv", index=False)


    # for group in hete_group_count_dict:
    #     print (group, hete_group_count_dict[group][0]/(sum(hete_group_count_dict[group])))
    
    data1 = []
    for gene in gene_list:
        for pop in hete_count_dict_pop:
            if gene in hete_count_dict_pop[pop]:
                if pop != 'Non-pop':
                    data1.append ([pop, gene, hete_count_dict_pop[pop][gene][0],hete_count_dict_pop[pop][gene][1],hete_count_dict_pop[pop][gene][0]/(sum(hete_count_dict_pop[pop][gene]))])
                else:
                    data1.append ([pop, gene,0, 0, 0])

    df = pd.DataFrame(data1, columns=['Pop','Gene', 'hete','hhomo','hete_freq'])
    df.to_csv("vdj_pop_hete_freq.csv", index=False)

def vdj_compare_hete_num(focus):
    df = pd.read_csv("vdj_snp_count.csv")

    # focus = 'AFR'
    df1 = df[df['Super Pop'] == focus]
    df2 = df[df['Super Pop'] != focus]

    # print (df1.shape, df2.shape)
    ## extract the gene list in df
    gene_list = set(df1['Gene'])
    data = []
    for gene in gene_list:
        hete1 = df1[df1['Gene'] == gene]['Hete_Variant_Num']
        hete2 = df2[df2['Gene'] == gene]['Hete_Variant_Num']
        ## compare hete1 and hete2 using t test
        t, p = scipy.stats.ttest_ind(hete1, hete2)
        ## compare fold change 
        if np.mean(hete2) == 0:
            fold_change = np.mean(hete1)
        else:
            fold_change = np.mean(hete1)/np.mean(hete2)

        # print (t, p)
        data.append([gene, t, p, np.mean(hete1), np.mean(hete2), len(hete1), len(hete2), len(hete1) + len(hete2), fold_change])
    # print (data)
    df = pd.DataFrame(data, columns=['Gene', 't', 'p_value', 'mean_hete1', 'mean_hete2', 'count1', 'count2', 'total', 'fold_change'])
    reject, pvals_corrected, _, alphacBonf = multipletests(list(df["p_value"]), alpha=0.05, method='bonferroni')
    df["p.adj"] = pvals_corrected
    df = df[df['p.adj'] < 0.01]
    ## sort df by p.adj
    df = df.sort_values(by=['p.adj'])
    print (df)
    print (df.shape)
    ## count how many has a fold change > 1, and the corresponding proportion
    enriched = df[df['fold_change'] > 1].shape[0]
    depleted = df[df['fold_change'] < 1].shape[0]

    print (focus, df[df['fold_change'] > 1].shape[0], df.shape[0], df[df['fold_change'] > 1].shape[0]/df.shape[0])
    df.to_csv(f"vdj_compare_hete_num_{focus}.csv", index=False)
    return enriched, depleted

def count_enriched_ratio():
    data = []
    for focus in ['AFR', 'EAS', 'EUR', 'SAS', 'AMR']:
        enriched, depleted = vdj_compare_hete_num(focus)
        data.append([focus, enriched, "Enriched"])
        data.append([focus, depleted, "Depleted"])
    df = pd.DataFrame(data, columns=['Pop', 'loci_num' ,'Type'])
    df.to_csv("vdj_compare_Enriched.csv", index=False)

def count_hete_freq():
    df = pd.read_csv("vdj_hete_freq.csv")
    ## add a new column, which is the sum of hete and homo
    df["total"] = df["Hete"] + df["Homo"]
    ## sort the df by hete_freq reversly
    df = df.sort_values(by=['Hete_freq'], ascending=False)
    print(df)

    
    


super_pop_file = "../hla/20131219.populations.tsv"
sample_pop_file = "../hla/20130606_sample_info.xlsx"
allele_file = "merged_samples.ig_tr.csv"

super_pop_dict = get_super_pop(super_pop_file)
sample_pop_dict = get_sample_pop(sample_pop_file)

# alleles_dict, alleles_gene_dict, alleles_sample_dict,pop_alleles_dict,sample_gene_dict \
#     = read_alleles(allele_file, super_pop_dict, sample_pop_dict, 10, 99)

# count_vdj(sample_gene_dict)
# count_enriched_ratio()
count_hete_freq()
