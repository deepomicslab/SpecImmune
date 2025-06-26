import pandas as pd
from collections import defaultdict
from scipy.stats import fisher_exact


def read_tab(csv, family="HLA"):
    sample_set = set()
    allele_sample_dict = defaultdict(dict)
    df = pd.read_csv(csv)
    for index, row in df.iterrows():
        sample = row['Sample']
        ## if row["One_guess"] is empty or not a string, skip
        if not isinstance(row["One_guess"], str) or row["One_guess"].strip()== "":
            continue
        # print (sample, row["One_guess"])
        ## get 1-field
        allele = row["One_guess"].split(":")[0]
        allele_sample_dict[allele][sample] = 1
        sample_set.add(sample)
    # print (len(allele_sample_dict), "alleles in", family)
    return allele_sample_dict, sample_set

def get_Contingency(allele_sample_dict, allele1, allele2, sample_set):
    """
    a: number of individuals with both alleles A and B

    b: A present, B absent

    c: A absent, B present

    d: neither present
    """
    ## Build 2x2 Contingency Tables
    a, b, c, d = 0, 0, 0, 0
    for sample in allele_sample_dict[allele1]:
        if sample in allele_sample_dict[allele2]:
            a += 1
        else:
            b += 1
    for sample in allele_sample_dict[allele2]:
        if sample not in allele_sample_dict[allele1]:
            c += 1
    d = len(sample_set) - a - b - c
    oddsratio, p_value = fisher_exact([[a, b], [c, d]])
    # print(a, b, c, d, oddsratio, p_value)
    return a, b, c, d, oddsratio, p_value

def association_test(HLA_allele_sample_dict, HLA_sample_set):
    allele_list = list(HLA_allele_sample_dict.keys())
    for i in range(len(allele_list)):
        for j in range(i + 1, len(allele_list)):
            allele1 = allele_list[i]
            allele2 = allele_list[j]
            if allele1 not in HLA_allele_sample_dict or allele2 not in HLA_allele_sample_dict:
                continue
            a, b, c, d, oddsratio, p_value = get_Contingency(HLA_allele_sample_dict, allele1, allele2, HLA_sample_set)
            if p_value < 0.01 and oddsratio > 1:
                print(f"{allele1} and {allele2}: a={a}, b={b}, c={c}, d={d}, oddsratio={oddsratio}, p_value={p_value}")

hla_csv = "/home/shuaiw/methylation/data/hla/genotypes/data/HLA_Table_S6.csv"
kir_csv = "/home/shuaiw/methylation/data/hla/genotypes/data/KIR_Table_S7.csv"
all_allele_sample_dict = {}
all_sample_set = set()
HLA_allele_sample_dict, HLA_sample_set = read_tab(hla_csv, family="HLA")
KIR_allele_sample_dict, KIR_sample_set = read_tab(kir_csv, family="KIR")
all_allele_sample_dict.update(HLA_allele_sample_dict)
all_allele_sample_dict.update(KIR_allele_sample_dict)
all_sample_set.update(HLA_sample_set)
all_sample_set.update(KIR_sample_set)
association_test(all_allele_sample_dict, all_sample_set)