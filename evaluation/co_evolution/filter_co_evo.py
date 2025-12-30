import pandas as pd    
from scipy.stats import binom
from statsmodels.stats.multitest import multipletests
import numpy as np

def replication_pvalue(x_obs, n_splits, p0):
    """
    Calculate binomial p-value for replication rate.
    
    Parameters
    ----------
    x_obs : int
        Number of observed replications
    n_splits : int
        Total number of splits tested (100)
    p0 : float
        Null probability of replication per split (0.05 for alpha=0.1)
    
    Returns
    -------
    float
        Upper-tail p-value: P(X >= x_obs | n, p0)
    """
    if x_obs == 0:
        return 1.0
    return binom.sf(x_obs - 1, n_splits, p0)

# Read results
total_df = pd.read_csv("./co_evo_pearson_results.csv")

# Filter for rows with replication data
total_df = total_df[total_df["split_sample_replication_rate"].notna()]
print(f"Total associations with replication data: {len(total_df)}")

# If replication_pvalue column doesn't exist, calculate it
if "replication_pvalue" not in total_df.columns:
    print("\nCalculating replication p-values...")
    # Assume 100 splits, alpha=0.1, so p0 = 0.5 * 0.1 = 0.05
    n_splits = 100
    alpha = 0.1
    p0 = 0.5 * alpha  # 0.05
    
    # Calculate number of replications and p-value
    total_df["replication_count"] = (total_df["split_sample_replication_rate"] * n_splits).round().astype(int)
    total_df["replication_pvalue"] = total_df["replication_count"].apply(lambda x: replication_pvalue(x, n_splits, p0))
    
    # Apply FDR correction
    _, pvals_corrected, _, _ = multipletests(total_df["replication_pvalue"], alpha=0.05, method='fdr_bh')
    total_df["replication_pvalue_FDR"] = pvals_corrected

# Filter by replication FDR
total_df_fdr = total_df[total_df["replication_pvalue_FDR"] < 0.05]
print(f"Associations with replication FDR < 0.05: {len(total_df_fdr)}")

# Also apply bootstrap filter
total_df_fdr = total_df_fdr[total_df_fdr["bootstrap_pval"] < 0.05]
print(f"After bootstrap filter (p < 0.05): {len(total_df_fdr)}")

# Save filtered results
total_df_fdr.to_csv("./co_evo_pearson_results.filtered.csv", index=False)

print(f"\nSummary statistics:")
print(f"  Mean replication rate: {total_df_fdr['split_sample_replication_rate'].mean():.3f}")
print(f"  Median replication rate: {total_df_fdr['split_sample_replication_rate'].median():.3f}")
print(f"  Min replication rate: {total_df_fdr['split_sample_replication_rate'].min():.3f}")
print(f"  Max replication rate: {total_df_fdr['split_sample_replication_rate'].max():.3f}")

    
