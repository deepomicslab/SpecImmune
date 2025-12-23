import pandas as pd
from collections import defaultdict
from scipy.stats import fisher_exact
from sklearn.metrics import mutual_info_score
import numpy as np
import re
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
## use seaborn for better aesthetics
import seaborn as sns
from scipy.stats import pearsonr
import gzip
import random
from cyvcf2 import VCF
from collections import defaultdict
from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression



def read_tab(csv, sample_pop_dict, family="HLA"):
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
        if row["Dataset"] != "1KGP":
            continue
        pop = sample_pop_dict[sample]
        
        gene = row["One_guess"].split("*")[0]
        ## get 1-field HLA
        if family == "HLA":
            field = row["One_guess"].split(":")
            allele = field[0]
            # if gene not in ["HLA-A", "HLA-B", "HLA-C", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "HLA-DRB1"]:
            #     continue
        if family == "KIR":
            field = row["One_guess"].split("*")
            allele = field[0] + "*" + field[1][:3]

        if pop not in allele_sample_dict[allele]:
            allele_sample_dict[allele][pop] = 0
        # print (pop)
        allele_sample_dict[allele][pop] += 1
        sample_set.add(sample)
        # if re.search("HLA-A\*01", allele):
        #     print (allele_sample_dict["HLA-A*01"])
    
    pop_num_dict = defaultdict(int)
    for sample in sample_set:
        pop = sample_pop_dict[sample]
        pop_num_dict[pop] += 1
    pop_list = sorted(list(pop_num_dict.keys()))
    print (pop_list)
    ## get frequency of each allele in each population
    allele_pop_freq_dict = defaultdict(dict)
    for allele in allele_sample_dict:
        allele_pop_freq_dict[allele] = {}
        for pop in pop_list:
            if pop not in allele_sample_dict[allele]:
                allele_pop_freq_dict[allele][pop] = 0
            else:
                allele_pop_freq_dict[allele][pop] = allele_sample_dict[allele][pop] / (pop_num_dict[pop]*2)

    return allele_pop_freq_dict, pop_list, allele_sample_dict, pop_num_dict

def read_cyp_csv(csv, sample_pop_dict):
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
        if row["Dataset"] != "1KGP":
            continue
        pop = sample_pop_dict[sample]
        
        gene = row["Locus"]
        allele = row["One_guess"].split(".")[0]

        if pop not in allele_sample_dict[allele]:
            allele_sample_dict[allele][pop] = 0
        allele_sample_dict[allele][pop] += 1
        sample_set.add(sample)
    
    pop_num_dict = defaultdict(int)
    for sample in sample_set:
        pop = sample_pop_dict[sample]
        pop_num_dict[pop] += 1
    pop_list = sorted(list(pop_num_dict.keys()))
    print (pop_list)
    ## get frequency of each allele in each population
    allele_pop_freq_dict = defaultdict(dict)
    for allele in allele_sample_dict:
        allele_pop_freq_dict[allele] = {}
        for pop in pop_list:
            if pop not in allele_sample_dict[allele]:
                allele_pop_freq_dict[allele][pop] = 0
            else:
                allele_pop_freq_dict[allele][pop] = allele_sample_dict[allele][pop] / (pop_num_dict[pop]*2)

    return allele_pop_freq_dict, pop_list, allele_sample_dict, pop_num_dict

def read_ig_tcr_csv(csv, sample_pop_dict, family="TR"):
    sample_set = set()
    allele_sample_dict = defaultdict(dict)
    
    df = pd.read_csv(csv)
    for index, row in df.iterrows():
        sample = row['Sample']
        if sample not in sample_pop_dict:
            continue
        if row["Dataset"] != "1KG":
            continue
        pop = sample_pop_dict[sample]
        
        gene = row["gene"]
        if gene[:2] != family:
            continue
        ## skip if allele is nan
        if pd.isna(row["allele_1"]) or pd.isna(row["allele_2"]):
            continue
        allele_list = [row["allele_1"], row["allele_2"]]
        for allele in allele_list:
            if pop not in allele_sample_dict[allele]:
                allele_sample_dict[allele][pop] = 0
            allele_sample_dict[allele][pop] += 1
        sample_set.add(sample)
    
    pop_num_dict = defaultdict(int)
    for sample in sample_set:
        pop = sample_pop_dict[sample]
        pop_num_dict[pop] += 1
    pop_list = sorted(list(pop_num_dict.keys()))
    print (pop_list)
    ## get frequency of each allele in each population
    allele_pop_freq_dict = defaultdict(dict)
    for allele in allele_sample_dict:
        allele_pop_freq_dict[allele] = {}
        for pop in pop_list:
            if pop not in allele_sample_dict[allele]:
                allele_pop_freq_dict[allele][pop] = 0
            else:
                allele_pop_freq_dict[allele][pop] = allele_sample_dict[allele][pop] / (pop_num_dict[pop]*2)

    return allele_pop_freq_dict, pop_list, allele_sample_dict, pop_num_dict

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

def plot(final_freq_dict, allele_1, allele_2, pop_list, super_pop_dict, r, p_adjusted):
    data2 = []
    for i in range(len(final_freq_dict[allele_1])):
        data2.append((final_freq_dict[allele_1][i], final_freq_dict[allele_2][i], pop_list[i], super_pop_dict[pop_list[i]]))
    df2 = pd.DataFrame(data2, columns=[f"{allele_1}_freq", f"{allele_2}_freq", "Population", "Super_Population"])
    ## plot the data

    sns.scatterplot(data=df2, x=f"{allele_1}_freq", y=f"{allele_2}_freq", hue="Super_Population")
    sns.regplot(data=df2, x=f"{allele_1}_freq", y=f"{allele_2}_freq", scatter=False, color='black')
    plt.xlabel(f"{allele_1} Frequency")
    plt.ylabel(f"{allele_2} Frequency")
    plt.title(f"Correlation: {r:.2f}, p.adj: {p_adjusted:.2e}")
    ## save the plot
    plt.savefig(f"plots/corr_{allele_1}_{allele_2}.pdf")
    ## clean the plot
    plt.clf()


def split_sample_replication(allele_sample_dict, pop_num_dict, pop_list, allele_pairs, n_splits=100, alpha=0.05):
    """
    Split-sample replication: partition individuals within each population,
    recompute frequencies, test edges in subset A, replicate in subset B.
    
    Returns replication rate for each allele pair.
    """
    replication_results = {pair: [] for pair in allele_pairs}
    
    for split_idx in range(n_splits):
        # Split samples within each population
        subset_A_dict = defaultdict(lambda: defaultdict(int))
        subset_B_dict = defaultdict(lambda: defaultdict(int))
        subset_A_count = defaultdict(int)
        subset_B_count = defaultdict(int)
        
        for allele, pop_counts in allele_sample_dict.items():
            for pop, count in pop_counts.items():
                # Randomly split count into two subsets
                count_A = np.random.binomial(count, 0.5)
                count_B = count - count_A
                subset_A_dict[allele][pop] = count_A
                subset_B_dict[allele][pop] = count_B
        
        # Recompute population sizes for each subset
        for pop in pop_list:
            pop_A_size = int(pop_num_dict[pop] * 0.5)
            pop_B_size = pop_num_dict[pop] - pop_A_size
            subset_A_count[pop] = pop_A_size
            subset_B_count[pop] = pop_B_size
        
        # Compute frequencies for subset A and B
        def compute_freqs(sample_dict, pop_count):
            freq_dict = {}
            for allele, pop_counts in sample_dict.items():
                freq_dict[allele] = {}
                for pop in pop_list:
                    count = pop_counts.get(pop, 0)
                    pop_size = pop_count.get(pop, 1)
                    freq_dict[allele][pop] = count / (pop_size * 2) if pop_size > 0 else 0
            return freq_dict
        
        freq_A = compute_freqs(subset_A_dict, subset_A_count)
        freq_B = compute_freqs(subset_B_dict, subset_B_count)
        
        # Test each allele pair
        for (allele1, allele2) in allele_pairs:
            if allele1 not in freq_A or allele2 not in freq_A:
                continue
            if allele1 not in freq_B or allele2 not in freq_B:
                continue
                
            freqs1_A = [freq_A[allele1][pop] for pop in pop_list]
            freqs2_A = [freq_A[allele2][pop] for pop in pop_list]
            freqs1_B = [freq_B[allele1][pop] for pop in pop_list]
            freqs2_B = [freq_B[allele2][pop] for pop in pop_list]
            
            # Skip if any zero frequencies
            if any(f == 0 for f in freqs1_A + freqs2_A + freqs1_B + freqs2_B):
                continue
            
            # Test in subset A
            corr_A, pval_A = pearsonr(freqs1_A, freqs2_A)
            
            # If significant in A, test in B
            if pval_A < alpha:
                corr_B, pval_B = pearsonr(freqs1_B, freqs2_B)
                # Concordant direction and nominal significance
                if np.sign(corr_A) == np.sign(corr_B) and pval_B < alpha:
                    replication_results[(allele1, allele2)].append(1)
                else:
                    replication_results[(allele1, allele2)].append(0)
    
    # Compute replication rate for each pair
    replication_rates = {}
    for pair, results in replication_results.items():
        if len(results) > 0:
            replication_rates[pair] = np.mean(results)
        else:
            replication_rates[pair] = np.nan
    
    return replication_rates


def bootstrap_correlation(freq1, freq2, n_bootstrap=1000):
    """
    Bootstrap resampling to assess correlation stability.
    Sample populations with replacement, keeping all allele frequencies within each population.
    """
    bootstrap_corrs = []
    n_pops = len(freq1)
    
    for _ in range(n_bootstrap):
        # Sample populations with replacement
        boot_indices = np.random.choice(n_pops, size=n_pops, replace=True)
        boot_freq1 = [freq1[i] for i in boot_indices]
        boot_freq2 = [freq2[i] for i in boot_indices]
        
        # Compute correlation on bootstrap sample
        corr, _ = pearsonr(boot_freq1, boot_freq2)
        bootstrap_corrs.append(corr)
    
    bootstrap_corrs = np.array(bootstrap_corrs)
    # Confidence intervals (95%)
    ci_lower = np.percentile(bootstrap_corrs, 2.5)
    ci_upper = np.percentile(bootstrap_corrs, 97.5)
    # Sign stability (proportion with same sign as observed)
    obs_corr = pearsonr(freq1, freq2)[0]
    sign_stability = np.mean(np.sign(bootstrap_corrs) == np.sign(obs_corr))
    # Bootstrap p-value (proportion crossing zero)
    bootstrap_pval = np.mean(np.sign(bootstrap_corrs) != np.sign(obs_corr))
    
    return ci_lower, ci_upper, sign_stability, bootstrap_pval


def permutation_test(HLA_allele_pop_freq_dict, KIR_allele_pop_freq_dict, pop_list, num_iterations=10000):
    """
    Perform a permutation test to calculate empirical p-values for the correlation between HLA and KIR allele frequencies.
    """
    permutation_correlations = []
    HLA_allele_list = list(HLA_allele_pop_freq_dict.keys())
    KIR_allele_list = list(KIR_allele_pop_freq_dict.keys())
    for i in range(num_iterations):
        ## randomly select one HLA allele and one KIR allele
        HLA_allele = random.sample(HLA_allele_list, 1)[0]
        KIR_allele = random.sample(KIR_allele_list, 1)[0]
        hla_freqs = [HLA_allele_pop_freq_dict[HLA_allele][pop] for pop in pop_list]
        kir_freqs = [KIR_allele_pop_freq_dict[KIR_allele][pop] for pop in pop_list]
        ## the freq in each pop should be >0
        if any(freq == 0 for freq in hla_freqs) or any(freq == 0 for freq in kir_freqs):
            continue
        ## get pearson correlation
        corr, _ = pearsonr(hla_freqs, kir_freqs)
        permutation_correlations.append(abs(corr))

    permutation_correlations = sorted(permutation_correlations)
    return permutation_correlations


def correlation_analysis(HLA_allele_pop_freq_dict, KIR_allele_pop_freq_dict, pop_list, corr_tag, empirical_correlations, pcs=None, 
                        HLA_sample_dict=None, KIR_sample_dict=None, pop_num_dict=None):
    corr_num = 0
    data = []
    final_freq_dict = {}
    allele_pairs_for_replication = []
    
    for HLA_allele in HLA_allele_pop_freq_dict:
        for KIR_allele in KIR_allele_pop_freq_dict:
            hla_freqs = [HLA_allele_pop_freq_dict[HLA_allele][pop] for pop in pop_list]
            kir_freqs = [KIR_allele_pop_freq_dict[KIR_allele][pop] for pop in pop_list]
            ## the freq in each pop should be >0
            if any(freq ==0 for freq in hla_freqs) or any(freq == 0 for freq in kir_freqs):
                continue
            final_freq_dict[HLA_allele] = hla_freqs
            final_freq_dict[KIR_allele] = kir_freqs
            ## get pearson correlation
            corr, pval = pearsonr(hla_freqs, kir_freqs)
            
            ## PC-adjusted p-value
            if pcs is not None:
                hla_arr = np.array(hla_freqs).reshape(-1, 1)
                kir_arr = np.array(kir_freqs).reshape(-1, 1)
                reg1 = LinearRegression().fit(pcs, hla_arr)
                reg2 = LinearRegression().fit(pcs, kir_arr)
                resid1 = hla_arr - reg1.predict(pcs)
                resid2 = kir_arr - reg2.predict(pcs)
                _, pval_pc = pearsonr(resid1.flatten(), resid2.flatten())
            else:
                pval_pc = pval
            
            if pval > 0.05 or pval_pc > 0.05:
                continue

            ## Bootstrap analysis
            ci_lower, ci_upper, sign_stability, bootstrap_pval = bootstrap_correlation(hla_freqs, kir_freqs, n_bootstrap=1000)
            
            ## calculate empirical p-value
            empirical_pval = sum(abs(corr) < abs(emp_corr) for emp_corr in empirical_correlations) / len(empirical_correlations)
            data.append((HLA_allele, KIR_allele, corr, pval, pval_pc, 
                        ci_lower, ci_upper, sign_stability, bootstrap_pval, 
                        corr_tag, empirical_pval))
            allele_pairs_for_replication.append((HLA_allele, KIR_allele))

    df = pd.DataFrame(data, columns=["HLA_allele", "KIR_allele", "Pearson_r", "p_value", "p_value_PC_adjusted",
                                      "bootstrap_CI_lower", "bootstrap_CI_upper", "bootstrap_sign_stability", "bootstrap_pval",
                                      "corr_tag", "empirical_p_value"])
    ## correct p-values

    p_values = df["p_value"].values
    rejected, pvals_corrected, _, _ = multipletests(p_values, alpha=0.05, method='fdr_bh')  ##multipletests(p_values, method='bonferroni')
    df["p_value_corrected"] = pvals_corrected
    ## keep only significant associations
    df = df[df["p_value_corrected"] < 0.05]
    print(f"Number of significant associations after correction: {len(df)}")
    ## sort by p_value_corrected
    df = df.sort_values(by="p_value_corrected", ascending=True)
    permutation_correlations = permutation_test(HLA_allele_pop_freq_dict, KIR_allele_pop_freq_dict, pop_list, 10000)
    ### add permutation p-value
    for index, row in df.iterrows():
        permutation_p_val = sum(abs(row["Pearson_r"]) < abs(emp_corr) for emp_corr in permutation_correlations) / len(permutation_correlations)
        df.at[index, "permutation_p_value"] = permutation_p_val
    df = df[df["permutation_p_value"] < 0.05]
    print (f"Number of significant associations after permutation test: {len(df)}")
    # print (df)
    df = df[df["empirical_p_value"] < 0.05]
    print (f"Number of significant associations after empirical test: {len(df)}")
    df = df[df["p_value_PC_adjusted"] < 0.05]
    print (f"Number of significant associations after PC-adjusted p-value: {len(df)}")
    
    # Split-sample replication
    if HLA_sample_dict is not None and KIR_sample_dict is not None and pop_num_dict is not None and len(df) > 0:
        print("\nPerforming split-sample replication (100 iterations)...")
        sig_pairs = [(row["HLA_allele"], row["KIR_allele"]) for _, row in df.iterrows()]
        combined_sample_dict = {**HLA_sample_dict, **KIR_sample_dict}
        replication_rates = split_sample_replication(combined_sample_dict, pop_num_dict, pop_list, sig_pairs, n_splits=100)
        df["split_sample_replication_rate"] = [replication_rates.get((row["HLA_allele"], row["KIR_allele"]), np.nan) 
                                                 for _, row in df.iterrows()]
        print(f"Mean replication rate: {df['split_sample_replication_rate'].mean():.3f}")
    
    print (df)
    # plot(final_freq_dict, "IGKV2-29*01", "HLA-DPA1*01", pop_list, super_pop_dict)
    return df, final_freq_dict

def read_snp_freq(chr1_vcf, sample_pop_dict, pop_list, gene_region_dict):
    """
    Parse a phased VCF file, randomly select 1000 SNPs,
    and compute ref allele frequencies per population.
    """
    vcf = VCF(chr1_vcf)
    samples = vcf.samples
    sample_to_index = {s: i for i, s in enumerate(samples)}

    # Filter samples for relevant populations
    pop_samples = {pop: [] for pop in pop_list}
    for sample in samples:
        pop = sample_pop_dict.get(sample)
        if pop in pop_list:
            pop_samples[pop].append(sample)

    print("Samples per population:")
    for pop, s_list in pop_samples.items():
        print(f"{pop}: {len(s_list)}")

    # Re-open VCF to reset iterator
    vcf = VCF(chr1_vcf)

    allele_pop_freq_dict = dict()

    for variant in vcf:
        ## randomly select with a rate of 1%
        if random.random() > 1:
            continue
        if not (variant.REF != '.' and variant.ALT[0] != '.'):
            continue
        ## only keep snp
        if len(variant.ALT) != 1 or len(variant.REF) != 1 or len(variant.ALT[0]) != 1:
            continue
        # print (variant.ID, variant.CHROM, variant.POS, variant.REF, variant.ALT)
        ### if vatiant is not in the gene region, skip
        if variant.CHROM not in gene_region_dict:
            continue
        in_gene_region = False
        for start, end in gene_region_dict[variant.CHROM]:
            if start <= variant.POS <= end:
                in_gene_region = True
                break
        if not in_gene_region:
            # print ("skip variant not in gene region", variant.ID, variant.CHROM, variant.POS)
            continue
        # if len(allele_pop_freq_dict) % 1000 == 0:
        #     print (af, len(allele_pop_freq_dict))
        hete_individuals = variant.INFO.get("AC_Het")
        individual_num = variant.INFO.get("AN") / 2  # Total individuals
        # print (hete_individuals,individual_num )
        hete_freq = hete_individuals / individual_num if individual_num > 0 else 0
        if hete_freq < 0.25:
            continue

        geno = variant.genotypes  # list of [a1, a2, phased]
        # build per-pop ref allele counts
        pop_ref_count = defaultdict(int)
        pop_total_alleles = defaultdict(int)
        for pop in pop_list:
            for sample in pop_samples[pop]:
                idx = sample_to_index[sample]
                alleles = geno[idx][:2]  # Get genotype (0 or 1)
                if None in alleles:
                    continue  # missing genotype
                ref_alleles = alleles.count(0)
                pop_ref_count[pop] += ref_alleles
                pop_total_alleles[pop] += 2

        # calculate ref frequency
        allele_pop_freq_dict[variant.ID] = {
            pop: (pop_ref_count[pop] / pop_total_alleles[pop]) if pop_total_alleles[pop] > 0 else None
            for pop in pop_list
        }
        # for pop in pop_list:
        #     print (f"{variant.ID} {pop}: {round(allele_pop_freq_dict[variant.ID][pop], 2)}, {pop_ref_count[pop]}, {pop_total_alleles[pop]}")

        if len(allele_pop_freq_dict) > 100:
            break
    print (len(allele_pop_freq_dict), "variants")
    return allele_pop_freq_dict

def empirical_test(chr1_allele_pop_freq_dict, chr10_allele_pop_freq_dict, pop_list, num_iterations=1000):
    ### randomly select 1000 alleles from chr1 and chr10, and calculate the correlation
    chr1_variants = list(chr1_allele_pop_freq_dict.keys())
    chr10_variants = list(chr10_allele_pop_freq_dict.keys())
    random.seed(2)  # for reproducibility
    selected_chr1_variants = random.sample(chr1_variants, num_iterations)
    selected_chr10_variants = random.sample(chr10_variants, num_iterations)
    empirical_correlations = []
    for i in range(num_iterations):
        chr1_freqs = [chr1_allele_pop_freq_dict[selected_chr1_variants[i]][pop] for pop in pop_list]
        chr10_freqs = [chr10_allele_pop_freq_dict[selected_chr10_variants[i]][pop] for pop in pop_list]
        if any(freq is None for freq in chr1_freqs) or any(freq is None for freq in chr10_freqs):
            continue
        corr, _ = pearsonr(chr1_freqs, chr10_freqs)
        empirical_correlations.append(abs(corr))
        if abs(corr) > 0.9:
            print (chr1_freqs)
            print (chr10_freqs)
            print (f"High correlation: {corr} {_} for {selected_chr1_variants[i]} and {selected_chr10_variants[i]}")
            print ("############")
            # ## plot the frequencies
            # plt.scatter(chr1_freqs, chr10_freqs)
            # plt.xlabel(f"Chr1 Frequencies ({selected_chr1_variants[i]})")
            # plt.ylabel(f"Chr10 Frequencies ({selected_chr10_variants[i]})")
            # plt.title(f"Correlation: {corr:.2f}")
            # plt.savefig(f"./empirical_corr_{selected_chr1_variants[i]}_{selected_chr10_variants[i]}.png")
            # plt.close() 
            # break
    empirical_correlations = sorted(empirical_correlations)
    ## plot box plot for the correlations
    plt.figure(figsize=(10, 6))
    sns.boxplot(x=empirical_correlations)
    plt.title("Empirical Correlations Distribution")
    plt.xlabel("Correlation Coefficient")
    plt.savefig("./empirical_correlations_boxplot.png")
    plt.close()
    print (empirical_correlations)
    return empirical_correlations

def read_gene_region(gtf):
    gene_region_dict = defaultdict(list)
    with open(gtf, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if fields[2] != "gene":
                continue
            chrom = fields[0]
            if chrom not in ["chr1", "chr10"]:
                continue
            if chrom not in gene_region_dict:
                gene_region_dict[chrom] = []
            start = min(int(fields[3]), int(fields[4]))
            end = max(int(fields[3]), int(fields[4]))
            gene_region_dict[chrom].append((start, end))

    return gene_region_dict

def alfred_empirical(alfred, num_iterations=10000):
    df = pd.read_csv(alfred, sep='\t')
    alfred_empirical_lists = df["corr_ref"].tolist()
    # alfred_empirical_lists = df["corr_ref"].tolist()
    if len(alfred_empirical_lists) < num_iterations:
        print(f"Warning: Not enough data points in ALFRED file. Using all {len(alfred_empirical_lists)} points.")
        num_iterations = len(alfred_empirical_lists)
    alfred_empirical_lists = random.sample(alfred_empirical_lists, num_iterations)
    return sorted(alfred_empirical_lists)

def output_raw_data(family_dict, pop_list):
    data = []
    for family, allele_pop_freq_dict in family_dict.items():
        for allele in allele_pop_freq_dict:
            allele_freqs = [allele_pop_freq_dict[allele][pop] for pop in pop_list]
            data.append((family, allele, *allele_freqs))
    df = pd.DataFrame(data, columns=["Family", "Allele"] + pop_list)
    df.to_csv(f"./allele_pop_freqs.tsv", sep='\t', index=False)


if __name__ == "__main__":

    hla_csv = "/home/shuaiw/methylation/data/hla/genotypes/data/HLA_Table_S6.csv"
    kir_csv = "/home/shuaiw/methylation/data/hla/genotypes/data/KIR_Table_S7.csv"
    cyp_csv = "/home/shuaiw/methylation/data/hla/genotypes/data/CYP_Table_S10.csv"
    ig_tcr_csv = "/home/shuaiw/methylation/data/hla/genotypes/data/IG_TR_Table_S9.csv"
    super_pop_file = "../1KG_ONT/hla/20131219.populations.tsv"
    sample_pop_file = "../1KG_ONT/hla/20130606_sample_info.xlsx"
    chr1_vcf = "/home/shuaiw/methylation/data/hla/phased_snps/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr1.filtered.shapeit2-duohmm-phased.vcf.gz"
    chr10_vcf = "/home/shuaiw/methylation/data/hla/phased_snps/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr10.filtered.shapeit2-duohmm-phased.vcf.gz"
    gtf = "/home/shuaiw/methylation/data/hla/gencode.v44.annotation.gtf"
    alfred = "ALFRED/pearson/alfred_200_random_interchrom_pairwise_population_corr.tsv"
    dbsnp = "dbSNP/dbSNP_pairwise_cor.tsv"
    gnomAD = "gnomAD/gene_pair_variant_pearson.tsv"
    gnomAD = "gnomAD/interchromosomal_gene_pairs_correlation.tsv"
    colors = ["#d9e6eb", "#9fc3d5", "#8f96bd", "#2a347a", "#d6d69b"]

    super_pop_dict = get_super_pop(super_pop_file)
    sample_pop_dict, pop_set, sample_super_pop_dict, super_pop_set = get_sample_pop(sample_pop_file, super_pop_dict)
    allele_count_dict = defaultdict(int)
    gene_region_dict = read_gene_region(gtf)
    print ("len(gene_region_dict)", len(gene_region_dict))

    pop_list = ['ACB', 'ASW', 'BEB', 'CDX', 'CEU', 'CHB', 'CHS', 'CLM', 'ESN', 'FIN', 'GBR', 'GIH', 'GWD', 'IBS', 'ITU', 'JPT', 'KHV', 'LWK', 'MSL', 'MXL', 'PEL', 'PJL', 'PUR', 'STU', 'TSI', 'YRI']
    # chr1_allele_pop_freq_dict = read_snp_freq(chr1_vcf, sample_pop_dict, pop_list, gene_region_dict)
    # chr10_allele_pop_freq_dict = read_snp_freq(chr10_vcf, sample_pop_dict, pop_list, gene_region_dict)
    # empirical_correlations = empirical_test(chr1_allele_pop_freq_dict, chr10_allele_pop_freq_dict, pop_list, 100)
    empirical_correlations = alfred_empirical(alfred, num_iterations=10000)

    TCR_pop_freq_dict, pop_list, TCR_sample_dict, TCR_pop_num_dict = read_ig_tcr_csv(ig_tcr_csv, sample_pop_dict, family="TR")
    IG_pop_freq_dict, pop_list, IG_sample_dict, IG_pop_num_dict = read_ig_tcr_csv(ig_tcr_csv, sample_pop_dict, family="IG")
    CYP_allele_pop_freq_dict, pop_list, CYP_sample_dict, CYP_pop_num_dict = read_cyp_csv(cyp_csv, sample_pop_dict)
    HLA_allele_pop_freq_dict, pop_list, HLA_sample_dict, HLA_pop_num_dict = read_tab(hla_csv, sample_pop_dict, family="HLA")
    KIR_allele_pop_freq_dict, pop_list, KIR_sample_dict, KIR_pop_num_dict = read_tab(kir_csv, sample_pop_dict, family="KIR")

    family_list = ["HLA", "KIR", "CYP", "IG", "TCR"]
    # family_list = ["CYP", "IG"]
    family_dict = {
        "HLA": HLA_allele_pop_freq_dict,
        "KIR": KIR_allele_pop_freq_dict,
        "CYP": CYP_allele_pop_freq_dict,
        "IG": IG_pop_freq_dict,
        "TCR": TCR_pop_freq_dict
    }

    # output_raw_data(family_dict, pop_list)

    # Compute PCs for PC-adjusted correlations
    all_freqs = []
    for family_name, allele_pop_freq_dict in family_dict.items():
        for allele, pop_freqs in allele_pop_freq_dict.items():
            freq_vector = [pop_freqs[pop] for pop in pop_list]
            if any(f > 0 for f in freq_vector):
                all_freqs.append(freq_vector)
    freq_matrix = np.array(all_freqs).T
    pca = PCA(n_components=min(5, freq_matrix.shape[0], freq_matrix.shape[1]))
    pcs = pca.fit_transform(freq_matrix)
    print(f"PCs computed: {pcs.shape}, variance explained: {pca.explained_variance_ratio_[:3]}")

    plot_info = []
    total_df = pd.DataFrame([])
    final_allele_freq_dict = {}
    for i in range(len(family_list)):
        for j in range(i+1, len(family_list)):
            family_1 = family_list[i]
            family_2 = family_list[j]
            corr_tag = f"{family_1}_vs_{family_2}"
            print(f"Correlation between {family_1} and {family_2}:")
            
            # Get sample dicts (for split-sample replication)
            sample_dicts = {"HLA": HLA_sample_dict, "KIR": KIR_sample_dict, "CYP": CYP_sample_dict, "IG": IG_sample_dict, "TCR": TCR_sample_dict}
            pop_num_dicts = {"HLA": HLA_pop_num_dict, "KIR": KIR_pop_num_dict, "CYP": CYP_pop_num_dict, "IG": IG_pop_num_dict, "TCR": TCR_pop_num_dict}
            sample_dict_1 = sample_dicts.get(family_1)
            sample_dict_2 = sample_dicts.get(family_2)
            pop_num = pop_num_dicts.get(family_1, HLA_pop_num_dict)
            
            df, final_freq_dict = correlation_analysis(
                family_dict[family_1], family_dict[family_2], pop_list, corr_tag, empirical_correlations, pcs=pcs,
                HLA_sample_dict=sample_dict_1, KIR_sample_dict=sample_dict_2, pop_num_dict=pop_num
            )
            ## get first row of df
            if not df.empty:
                print(df.iloc[0])
            plot_info.append([corr_tag, df.iloc[0]["HLA_allele"], df.iloc[0]["KIR_allele"], df.iloc[0]["Pearson_r"], df.iloc[0]["p_value_corrected"]])
            total_df = pd.concat([total_df, df], ignore_index=True)
            final_allele_freq_dict = {**final_allele_freq_dict, **final_freq_dict}
        #     if j > 1:
        #         break
        # break

            
    # save the total df to csv
    total_df.to_csv(f"./co_evo_pearson_results.csv", index=False)

    
    # # Create a 5x2 grid for subplots
    # fig, axes = plt.subplots(5, 2, figsize=(12, 24))
    # axes = axes.flatten()
    # super_pop_order = sorted(list(super_pop_set))


    # for idx, info in enumerate(plot_info):
    #     if idx >= 10:
    #         break
    #     corr_tag, allele_1, allele_2, r, p_adjusted = info
    #     print(f"Plotting {corr_tag} for {allele_1} and {allele_2} with r={r:.2f}, p.adj={p_adjusted:.2e}")
    #     # Prepare data
    #     data2 = []
    #     for i in range(len(final_allele_freq_dict[allele_1])):
    #         data2.append((final_allele_freq_dict[allele_1][i], final_allele_freq_dict[allele_2][i], pop_list[i], super_pop_dict[pop_list[i]]))
    #     df2 = pd.DataFrame(data2, columns=[f"{allele_1}_freq", f"{allele_2}_freq", "Population", "Super_Population"])
    #     ## sort the df2 by super_pop_order
    #     df2["Super_Population"] = pd.Categorical(df2["Super_Population"], categories=super_pop_order, ordered=True)
    #     df2 = df2.sort_values("Super_Population")
    #     ax = axes[idx]
    #     if idx == 0:
    #         scatter = sns.scatterplot(data=df2, x=f"{allele_1}_freq", y=f"{allele_2}_freq", hue="Super_Population", ax=ax, palette=colors)
    #     else:
    #         scatter = sns.scatterplot(data=df2, x=f"{allele_1}_freq", y=f"{allele_2}_freq", hue="Super_Population",ax=ax, legend=False, palette=colors)
    #     sns.regplot(data=df2, x=f"{allele_1}_freq", y=f"{allele_2}_freq", scatter=False, color='black', ax=ax)
    #     ax.set_xlabel(f"{allele_1} Frequency")
    #     ax.set_ylabel(f"{allele_2} Frequency")
    #     ax.set_title(f"Correlation: {r:.2f}, p.adj: {p_adjusted:.2e}")



    # plt.savefig("plots/corr_all_subplots.pdf")
    # plt.clf()



    # plot(final_allele_freq_dict, "HLA-DPB1*04", "KIR2DL4*008", pop_list, super_pop_dict)

    # correlation_analysis(IG_pop_freq_dict, HLA_allele_pop_freq_dict, pop_list)
    # correlation_analysis(KIR_allele_pop_freq_dict, CYP_allele_pop_freq_dict, pop_list)

