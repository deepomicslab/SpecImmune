import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, ttest_ind
import numpy as np

# Step 1: Define gene classes
ordered_loci = ['KIR2DL1', 'KIR2DL2', 'KIR2DL3', 'KIR2DL4', 'KIR2DL5A', 'KIR2DL5B', 'KIR2DP1', 'KIR2DS1', 'KIR2DS2', 'KIR2DS3', 'KIR2DS4', 'KIR2DS5', 'KIR3DL1', 'KIR3DL2', 'KIR3DL3', 'KIR3DP1', 'KIR3DS1']

# Step 2: Read the data from the file
file_path = 'superpopulation_diversity_df_KIR.csv'  # specify your file path
df = pd.read_csv(file_path)

# Step 3: Pivot the DataFrame to get the correct format for the heatmap
df_pivot = df.pivot(index='Superpopulation', columns='Locus', values='Diversity')

# Step 4: Reorder the columns (Loci) based on the gene classes
df_pivot = df_pivot[ordered_loci]

# Step 5: Calculate ranks for each gene (loci) across superpopulations (column-wise ranking)
df_ranks = df_pivot.rank(axis=0, method='min', ascending=False)

# Step 6: Separate EAS and non-EAS populations
eas_population = ['EAS']  # Modify this list with the actual East Asian population name(s)
non_eas_population = df_ranks.index.difference(eas_population)

# Extract the 1x17 rank vector for EAS
eas_ranks_vector = df_ranks.loc[eas_population].values.flatten()

# Extract the rank matrix for non-EAS
non_eas_ranks_matrix = df_ranks.loc[non_eas_population].values

# Flatten the non-EAS matrix into a 1D vector
non_eas_ranks_vector = non_eas_ranks_matrix.flatten()

# Step 7: Perform the Mann-Whitney U test or t-test to compare the EAS vector vs. non-EAS vector
# Mann-Whitney U test (non-parametric)
stat, mannwhitney_p_value = mannwhitneyu(eas_ranks_vector, non_eas_ranks_vector, alternative='two-sided')

# T-test (parametric)
t_stat, ttest_p_value = ttest_ind(eas_ranks_vector, non_eas_ranks_vector, equal_var=False)

# Print the results
print("\nOverall significance (Mann-Whitney U Test) p-value:", mannwhitney_p_value)
print("Overall significance (T-Test) p-value:", ttest_p_value)

# Step 8: Create a DataFrame for plotting
eas_labels = ['EAS'] * len(eas_ranks_vector)
non_eas_labels = ['Non-EAS'] * len(non_eas_ranks_vector)

ranks_df = pd.DataFrame({
    'Rank': np.concatenate([eas_ranks_vector, non_eas_ranks_vector]),
    'Population': eas_labels + non_eas_labels
})

# Step 9: Plot Box Plots using Seaborn and add T-Test p-value
plt.figure(figsize=(10, 6))

sns.boxplot(x='Population', y='Rank', data=ranks_df, palette='Set2', width=0.4)

# Adding the T-Test p-value on the plot in scientific notation
plt.title('Distribution of KIR locus diversity ranks for EAS vs Non-EAS populations', weight='bold')
plt.ylabel('Rank', weight='bold')
plt.xlabel('Population', weight='bold')
plt.xticks(weight='bold')

# Adjust the y-axis to make room for the p-value text
plt.ylim([0, ranks_df['Rank'].max() + 1])

# Invert the y-axis to show large values at the bottom and small values at the top
plt.gca().invert_yaxis()

# Position the p-value text on the plot in scientific notation
plt.text(0.5, ranks_df['Rank'].min() - 0.5, f'T-Test p-value: {ttest_p_value:.2e}',
         horizontalalignment='center', size='medium', color='black', weight='bold')

plt.show()