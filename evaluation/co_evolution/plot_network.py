import pandas as pd


corr_csv = "co_evo_pearson_results.csv"
df = pd.read_csv(corr_csv)
df = df[df["corr_tag"] == "HLA_vs_KIR"]


pair_dict = {}
for i, row in df.iterrows():
    gene_1 = row["HLA_allele"].split("*")[0]
    gene_2 = row["KIR_allele"].split("*")[0]
    r = row["Pearson_r"]

    tag = f"{gene_1}_{gene_2}"
    if tag not in pair_dict:
        pair_dict[tag] = []
    pair_dict[tag].append(r)

for tag, r_list in pair_dict.items():
    print (tag, r_list)