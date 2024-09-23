import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# Step 1: Define gene classe

# Combine all classes into a single ordered list for the x-axis
ordered_loci = ['KIR2DL1', 'KIR2DL2', 'KIR2DL3', 'KIR2DL4', 'KIR2DL5A', 'KIR2DL5B', 'KIR2DP1', 'KIR2DS1', 'KIR2DS2', 'KIR2DS3', 'KIR2DS4', 'KIR2DS5', 'KIR3DL1', 'KIR3DL2', 'KIR3DL3', 'KIR3DP1', 'KIR3DS1']


# Step 2: Read the data from the file
file_path = 'superpopulation_diversity_df_KIR.csv'  # specify your file path
df = pd.read_csv(file_path)

# Step 3: Pivot the DataFrame to get the correct format for the heatmap
df_pivot = df.pivot(index='Superpopulation', columns='Locus', values='Diversity')


# Remove "HLA-" from the Locus names
df_pivot.columns = df_pivot.columns.str.replace('HLA-', '', regex=False)

# Step 4: Reorder the columns (Loci) based on the gene classes
df_pivot = df_pivot[ordered_loci]


# Step 5: Calculate average diversity for each gene and each population
# Average diversity for each gene (column-wise mean)
gene_avg_diversity = df_pivot.mean(axis=0)

# Average diversity for each population (row-wise mean)
population_avg_diversity = df_pivot.mean(axis=1)

# Print the results
print("Average diversity for each gene (locus):")
print(gene_avg_diversity)

print("\nAverage diversity for each population:")
print(population_avg_diversity)

# Step 6: Create a custom colormap
custom_cmap = LinearSegmentedColormap.from_list("custom_cmap", ["lightgray", "#C6EDCF", "#22C9BF", "#49C7DA", "#006BB2", '#002786'])

# Step 7: Plot the heatmap
plt.figure(figsize=(20, 10))
sns.heatmap(df_pivot, annot=True, linewidths=1, annot_kws={"size": 12}, cbar_kws={}, cmap=custom_cmap)

# Customize the plot
plt.title('', fontsize=20, fontweight='bold')
plt.xlabel('Locus', fontsize=16, fontweight='bold')
plt.ylabel('Population', fontsize=16, fontweight='bold')

# Adjust the tick labels size
plt.xticks(fontsize=12, rotation=45, ha='right')
plt.yticks(fontsize=12)

# Save the plot as svg
plt.savefig('diversity-heatmap.svg', format='svg', bbox_inches='tight', dpi=600)

# Show the plot
plt.show()
