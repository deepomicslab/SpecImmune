import pandas as pd

import sys, os
sys.path.insert(0, sys.path[0]+'/../')
from four_field_compare import convert_field


## read the file using pandas, sep is \t
def get_super_pop(super_pop_file):
    df = pd.read_csv(super_pop_file, sep='\t')
    super_pop_dict = {}
    for index, row in df.iterrows():
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

def read_alleles(allele_file, super_pop_dict, sample_pop_dict):
    df = pd.read_csv(allele_file)
    alleles_dict = {}
    for index, row in df.iterrows():
        if row['Sample'] not in sample_pop_dict:
            # print ("no pop info for sample", row['Sample'])
            continue
        # else:
        #     print (row['Sample'], sample_pop_dict[row['Sample']])
        pop = sample_pop_dict[row['Sample']]
        super_pop = super_pop_dict[pop]
        if super_pop not in alleles_dict:
            alleles_dict[super_pop] = {}
        if row['Locus'] not in alleles_dict[super_pop]:
            alleles_dict[super_pop][row['Locus']] = []
        ## check if empty
        if row['Genotype'] != "NA" and row['Genotype'] != "" and not pd.isna(row['Genotype']):
            genotype = row['Genotype']
            genotype = convert_field([genotype], 4)[0]
            alleles_dict[super_pop][row['Locus']].append(genotype)
    return alleles_dict

def count_freq(list):
    ## count the frequency of each allele
    freq_dict = {}
    total = 0
    for allele in list:
        if allele not in freq_dict:
            freq_dict[allele] = 0
        total += 1
        freq_dict[allele] += 1
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
        if super_pop not in freq_dict:
            freq_dict[super_pop] = {}
        all_alleles = []
        for locus in alleles_dict[super_pop]:
            all_alleles += alleles_dict[super_pop][locus]
        all_pop_alleles += all_alleles
        freq_dict = count_freq(all_alleles)
        super_pop_freq_dist[super_pop] = freq_dict

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


super_pop_file = "hla/20131219.populations.tsv"
sample_pop_file = "hla/20130606_sample_info.xlsx"
allele_file = "hla/speclong_res_merged_samples.csv"
freq_file = "hla/hla_freq.csv"
super_pop_dict = get_super_pop(super_pop_file)
sample_pop_dict = get_sample_pop(sample_pop_file)
alleles_dict = read_alleles(allele_file, super_pop_dict, sample_pop_dict)
count_alleles(alleles_dict, freq_file)