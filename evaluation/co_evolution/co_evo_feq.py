import pandas as pd
from collections import defaultdict
from scipy.stats import fisher_exact


def read_tab(csv, sample_pop_dict, pop, family="HLA"):
    sample_set = set()
    allele_sample_dict = defaultdict(dict)
    df = pd.read_csv(csv)
    for index, row in df.iterrows():
        sample = row['Sample']
        ## if row["One_guess"] is empty or not a string, skip
        if not isinstance(row["One_guess"], str) or row["One_guess"].strip()== "":
            continue
        if sample not in sample_pop_dict:
            continue
        if sample_pop_dict[sample] != pop:
            continue
        # print (sample, row["One_guess"])
        ## get 1-field HLA
        if family == "HLA":
            allele = row["One_guess"].split(":")[0]
        if family == "KIR":
            allele = row["One_guess"][:11]
        if allele.split("*")[0] in ["MICA", "MICB", "TAP2"]:
            allele = "HLA-" + allele
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

def association_test(all_allele_sample_dict, all_sample_set):
    allele_list = list(all_allele_sample_dict.keys())
    # allele_list = allele_list[:5]
    data = []
    for i in range(len(allele_list)):
        for j in range(i + 1, len(allele_list)):
            allele1 = allele_list[i]
            allele2 = allele_list[j]
            if allele1.split("*")[0] == allele2.split("*")[0]:  ## remove same gene
                continue
            
            if allele1[:3] == allele2[:3]:  ## remove same family
                continue
            if allele1 not in all_allele_sample_dict or allele2 not in all_allele_sample_dict:
                continue
            a, b, c, d, oddsratio, p_value = get_Contingency(all_allele_sample_dict, allele1, allele2, all_sample_set)
            data.append((allele1, allele2, a, b, c, d, oddsratio, p_value))
            # if p_value < 0.01 and oddsratio > 1:
            #     print(f"{allele1} and {allele2}: a={a}, b={b}, c={c}, d={d}, oddsratio={oddsratio}, p_value={p_value}")
    print (f"start correction..")
    ## perform multiple testing correction
    from statsmodels.stats.multitest import multipletests
    df = pd.DataFrame(data, columns=["Allele1", "Allele2", "a", "b", "c", "d", "OddsRatio", "p_value"])
    ## use rejected, pvals_corrected, _, _ = multipletests(p_values, alpha=0.05, method='fdr_bh')
    p_values = df["p_value"].values
    rejected, pvals_corrected, _, _ = multipletests(p_values, alpha=0.05, method='fdr_bh')

    ## add p_value_corrected to df
    df["p_value_corrected"] = pvals_corrected
    ## sort by p_value_corrected
    df = df.sort_values(by="p_value_corrected")
    ## print top 10 associations
    print(df.head(10))
    ## filter by p_value_corrected < 0.05
    df_filtered = df[df["p_value_corrected"] < 0.05]
    print(f"Number of associations with p_value_corrected < 0.05: {len(df_filtered)}")
    ## save to file


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
def get_sample_pop(sample_pop_file, super_pop_dict):
    pop_set = set()
    super_pop_set = set()
    # read excel file into df
    df = pd.read_excel(sample_pop_file)
    sample_pop_dict = {}
    sample_super_pop_dict = {}
    for index, row in df.iterrows():
        sample_pop_dict[row['Sample']] = row['Population']
        super_pop = super_pop_dict.get(row['Population'], None)
        pop_set.add(row['Population'])
        super_pop_set.add(super_pop)
        sample_super_pop_dict[row['Sample']] = super_pop
        # print (f"{row['Sample']} -> {row['Population']}")
    # print (sample_pop_dict)
    return sample_pop_dict, pop_set, sample_super_pop_dict, super_pop_set

hla_csv = "/home/shuaiw/methylation/data/hla/genotypes/data/HLA_Table_S6.csv"
kir_csv = "/home/shuaiw/methylation/data/hla/genotypes/data/KIR_Table_S7.csv"
super_pop_file = "../1KG_ONT/hla/20131219.populations.tsv"
sample_pop_file = "../1KG_ONT/hla/20130606_sample_info.xlsx"

super_pop_dict = get_super_pop(super_pop_file)
sample_pop_dict, pop_set, sample_super_pop_dict, super_pop_set = get_sample_pop(sample_pop_file, super_pop_dict)


for pop in super_pop_set:
    all_allele_sample_dict = {}
    all_sample_set = set()
    HLA_allele_sample_dict, HLA_sample_set = read_tab(hla_csv, sample_super_pop_dict, pop, family="HLA")
    KIR_allele_sample_dict, KIR_sample_set = read_tab(kir_csv, sample_super_pop_dict, pop, family="KIR")
    all_allele_sample_dict.update(HLA_allele_sample_dict)
    all_allele_sample_dict.update(KIR_allele_sample_dict)
    all_sample_set.update(HLA_sample_set)
    all_sample_set.update(KIR_sample_set)
    print (pop, "pop...", len(all_allele_sample_dict), "alleles,",
           len(all_sample_set), "samples")
    association_test(all_allele_sample_dict, all_sample_set)