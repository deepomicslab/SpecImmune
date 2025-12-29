import pandas as pd    
    
total_df = pd.read_csv(f"./co_evo_pearson_results.csv")
## filter total_df by split_sample_replication_rate exists and split_sample_replication_rate > 0.5
total_df = total_df[total_df["split_sample_replication_rate"].notna()]
total_df = total_df[total_df["split_sample_replication_rate"] >= 0.7]
total_df = total_df[total_df["bootstrap_pval"] < 0.05]
print (f"Number of significant associations after bootstrap analysis: {len(total_df)}")
total_df.to_csv(f"./co_evo_pearson_results.filtered.csv", index=False)

    
