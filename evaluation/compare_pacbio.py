import pandas as pd

def compare_digit_level(allele1, allele2, digit=8):
    allele_str = allele1.split("*")[1]
    truth_allele_str = allele2.split("*")[1]
    return truth_allele_str in allele_str

def has_intersection(alleles, truth_alleles, digit=8):
    return any(compare_digit_level(allele, truth_allele, digit) for allele in alleles for truth_allele in truth_alleles)

def process_hla_data(hlala_file, record_file, top_n=5):
    genes = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]
    samples = ["HG00514", "HG00731", "HG00732", "HG00733", "NA19238", "NA19239", "NA19240"]

    hlala_dict = {sample: {gene: [[], []] for gene in genes} for sample in samples}

    # Parse HLA-LA result file
    with open(hlala_file, "r") as hlala_f:
        for idx, line in enumerate(hlala_f):
            if idx == 0:  # Skip header
                continue
            parts = line.strip().split('\t')
            sample = parts[0].split("_")[0]
            
            if sample not in samples:
                continue

            gene_visited = {gene: 0 for gene in genes}
            for pid, part in enumerate(parts[1:], start=1):
                if part == "":
                    continue
                for allele in part.split(";"):
                    if allele == "":
                        continue
                    allele_gene = allele.split("*")[0]
                    if gene_visited[allele_gene] == 0:
                        hlala_dict[sample][allele_gene][0].append(allele)
                    else:
                        hlala_dict[sample][allele_gene][1].append(allele)
                    gene_visited[allele_gene] = 1

    # Parse truth data
    df = pd.read_csv(record_file, sep='\s+', header=None, names=['sample', 'gene', 'hap_idx', 'allele', 'length', 'identity'])

    gene_right_count_dict = {}

    for sample in samples:
        for gene in genes:
            all_identity_h1 = list(set(df.loc[(df['sample'] == sample) & (df['gene'] == gene) & (df['hap_idx'] == 'h1'), 'identity']))
            all_identity_h1.sort(reverse=True)
            all_identity_h1 = all_identity_h1[:top_n]
            
            truth_allele1 = set(df.loc[(df['sample'] == sample) & (df['gene'] == gene) & (df['hap_idx'] == 'h1') & (df['identity'].isin(all_identity_h1)), 'allele'])
            
            all_identity_h2 = list(set(df.loc[(df['sample'] == sample) & (df['gene'] == gene) & (df['hap_idx'] == 'h2'), 'identity']))
            all_identity_h2.sort(reverse=True)
            all_identity_h2 = all_identity_h2[:top_n]
            
            truth_allele2 = set(df.loc[(df['sample'] == sample) & (df['gene'] == gene) & (df['hap_idx'] == 'h2') & (df['identity'].isin(all_identity_h2)), 'allele'])
            
            hla_allele1 = hlala_dict[sample][gene][0]
            hla_allele2 = hlala_dict[sample][gene][1]
            
            inter_0_0_cnt = int(has_intersection(hla_allele1, truth_allele1))
            inter_0_1_cnt = int(has_intersection(hla_allele1, truth_allele2))
            inter_1_0_cnt = int(has_intersection(hla_allele2, truth_allele1))
            inter_1_1_cnt = int(has_intersection(hla_allele2, truth_allele2))
            
            max_inter_len = max(inter_0_0_cnt + inter_1_1_cnt, inter_1_0_cnt + inter_0_1_cnt)
            
            if max_inter_len < 2:
                print(sample, (truth_allele1, truth_allele2), hlala_dict[sample][gene], max_inter_len)
            
            gene_right_count_dict[(sample, gene)] = max_inter_len

    # Calculate accuracy for each gene
    gene_accuracy = {gene: sum(gene_right_count_dict[(sample, gene)] for sample in samples) / (2 * len(samples)) for gene in genes}

    # Print accuracy for each gene
    for gene, acc in gene_accuracy.items():
        print(f"{gene}: {acc}")
    print(gene_accuracy)

# Example usage
hlala_file = "HLA-LA.merge.result.txt"
record_file = "record_alignment_file.txt"
process_hla_data(hlala_file, record_file)