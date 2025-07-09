import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Read the results file
df = pd.read_csv("alfred_200_random_interchrom_pairwise_population_corr.tsv", sep='\t')

# Keep only valid correlation values
df = df[df['corr_ref'].notnull()]

# Take the absolute value of correlation coefficients
df['abs_corr_ref'] = df['corr_ref'].abs()

plt.figure(figsize=(5, 7))
sns.boxplot(y=df['abs_corr_ref'], color='#3486eb', width=0.3)

plt.ylabel('Absolute ref allele frequency correlation (|Pearson r|)')
plt.title('Distribution of Absolute Interchromosomal SNP Pairwise Correlations')
plt.grid(axis='y', linestyle='--', alpha=0.5)
plt.tight_layout()
plt.savefig('pairwise_abs_corr_boxplot.png', dpi=150)
plt.show()