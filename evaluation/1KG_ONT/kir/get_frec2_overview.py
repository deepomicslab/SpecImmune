import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import os

# ---------------------------- Configuration & Setup ----------------------------

# Ensure that the script is called with the required output directory argument
if len(sys.argv) < 2:
    print("Usage: python script.py <output_directory>")
    sys.exit(1)

outdir = sys.argv[1]
if not os.path.exists(outdir):
    os.makedirs(outdir)

# ---------------------------- Helper Functions ----------------------------

def get_population(df, sample_name):
    """
    Given a DataFrame and a sample name, return the population.

    Parameters:
    df (pd.DataFrame): The DataFrame containing the sample information.
    sample_name (str): The name of the sample to look up.

    Returns:
    str: The population of the sample.
    """
    sample_row = df[df['Sample'] == sample_name]
    if sample_row.empty:
        return "Unknown"
    return sample_row['Population'].values[0]

def calculate_shannon_diversity(frequencies):
    """
    Calculate the Shannon diversity index for a set of frequencies.

    Parameters:
    frequencies (pd.Series): A series of allele frequencies.

    Returns:
    float: Shannon diversity index.
    """
    frequencies = frequencies[frequencies > 0]
    return -np.sum(frequencies * np.log(frequencies))

# ---------------------------- Data Loading & Preprocessing ----------------------------

# Step 1: Read the gene read counts file into a DataFrame
try:
    gene_reads_df = pd.read_csv('./read_depth_kirs.csv', index_col=0)
except FileNotFoundError:
    print("Error: 'read_depth_kirs.csv' not found.")
    sys.exit(1)

print("Gene Reads DataFrame:")
print(gene_reads_df.head())

# Step 2: Set the reads cutoff values
read_cutoffs = [10]  # Add more cutoff values if needed

# Step 3: Read the population and region data from the Excel file
try:
    metadata_df = pd.read_excel('./20130606_sample_info.xlsx', engine='openpyxl')
except FileNotFoundError:
    print("Error: '20130606_sample_info.xlsx' not found.")
    sys.exit(1)

# Ensure consistent column names by replacing spaces with underscores
metadata_df.columns = metadata_df.columns.str.replace(' ', '_')

# Step 4: Read the allele/genotype data
try:
    allele_df = pd.read_csv('merged_samples.csv')
except FileNotFoundError:
    print("Error: 'merged_samples.csv' not found.")
    sys.exit(1)

# ---------------------------- Data Transformation ----------------------------

# Add population information to the allele DataFrame
allele_df['Population'] = allele_df['Sample'].apply(lambda x: get_population(metadata_df, x))

# Split alleles in the Genotype column by ';' and take the first allele
allele_df['Alleles'] = allele_df['Genotype'].str.split(';').str[0]

# Save the updated allele DataFrame to a CSV file
allele_df.to_csv('allele_df.csv', index=False)
print("Allele DataFrame saved to 'allele_df.csv'.")

# Group by Population, Locus (gene), and Allele to count occurrences
allele_counts = allele_df.groupby(['Population', 'Locus', 'Alleles']).size().reset_index(name='Count')

# Add a "Total" column which is the sum of counts for the same Population and Locus
allele_counts['Total'] = allele_counts.groupby(['Population', 'Locus'])['Count'].transform('sum')

# Calculate the frequency of each allele by dividing the count by the total
allele_counts['Frequency'] = allele_counts['Count'] / allele_counts['Total']

# Save the allele counts DataFrame to a CSV file
allele_counts.to_csv('allele_counts.csv', index=False)
print("Allele counts saved to 'allele_counts.csv'.")

# Calculate the Shannon diversity index for each Locus in each Population
diversity_df = allele_counts.groupby(['Population', 'Locus'])['Frequency'].apply(calculate_shannon_diversity).reset_index(name='Diversity')

# Save the diversity DataFrame to a CSV file
diversity_df.to_csv('diversity_df.csv', index=False)
print("Diversity DataFrame saved to 'diversity_df.csv'.")

# ---------------------------- Plotting Configuration ----------------------------

# Set the style for the plots
sns.set(style="white")

# Define a colormap
cmap = plt.cm.viridis

# ---------------------------- Plotting Per Population ----------------------------

# Get the list of unique populations
populations = allele_counts['Population'].unique()

# Debugging: Print the unique populations
print(f"Unique populations found (excluding 'Unknown'): {populations[populations != 'Unknown']}")

# ---------------------------- Plotting for All Samples ----------------------------

# Aggregate allele counts across all populations
allele_counts_all = allele_counts.groupby(['Locus', 'Alleles']).agg({'Count': 'sum'}).reset_index()

# Calculate total counts per locus
allele_counts_all['Total'] = allele_counts_all.groupby('Locus')['Count'].transform('sum')

# Calculate frequency
allele_counts_all['Frequency'] = allele_counts_all['Count'] / allele_counts_all['Total']

# Calculate Shannon diversity for all samples per locus
diversity_all = allele_counts_all.groupby('Locus')['Frequency'].apply(calculate_shannon_diversity).reset_index(name='Diversity')

# Save the aggregated allele counts and diversity data for all samples
allele_counts_all.to_csv('allele_counts_all.csv', index=False)
diversity_all.to_csv('diversity_all_df.csv', index=False)
print("Aggregated allele counts saved to 'allele_counts_all.csv'.")
print("Aggregated diversity data saved to 'diversity_all_df.csv'.")

# Pivot the data for plotting
pivot_df_all = allele_counts_all.pivot(index='Alleles', columns='Locus', values='Frequency').fillna(0)

# Sort alleles by their overall frequency for better visualization
pivot_df_all = pivot_df_all.reindex(pivot_df_all.sum(axis=1).sort_values(ascending=False).index)

# Debugging: Print the sorted alleles for all samples
print("\nTop alleles sorted by total frequency for all samples:")
print(pivot_df_all.sum(axis=1).head())

# Debugging: Ensure that diversity is sorted correctly
diversity_all_sorted = diversity_all.sort_values(by='Diversity', ascending=False)
sorted_loci_all = diversity_all_sorted['Locus'].tolist()
pivot_df_all = pivot_df_all[sorted_loci_all]

print(f"Order of loci for all samples plot (sorted by diversity): {sorted_loci_all}")

# Debugging: Print diversity values to ensure correct sorting
print("\nShannon diversity values for all samples (sorted):")
print(diversity_all_sorted)

# ---------------------------- Define Global Normalization Based on Frequencies ----------------------------

# **Step 1: Calculate global min and max frequencies across all alleles and loci**
global_min_frequency = allele_counts_all['Frequency'].min()
global_max_frequency = allele_counts_all['Frequency'].max()

print(f"Global Minimum Allele Frequency: {global_min_frequency}")
print(f"Global Maximum Allele Frequency: {global_max_frequency}")

# **Step 2: Define normalization based on global min and max frequencies**
norm = plt.Normalize(vmin=global_min_frequency, vmax=global_max_frequency)

# ---------------------------- Plotting Shannon Diversity and Allele Frequencies ----------------------------

# Plot Shannon diversity alongside allele frequencies for all samples
fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(120, 90), gridspec_kw={'height_ratios': [1, 10]})

# --- Top plot: Shannon diversity line plot ---
diversity_values_all = diversity_all_sorted['Diversity'].values
ax1.plot(sorted_loci_all, diversity_values_all, marker='o', color='red', linewidth=6, markersize=30)

# Add title and labels
ax1.set_title('Shannon diversity across KIR genes', fontsize=96, fontweight='bold')
ax1.set_ylabel('Diversity', fontsize=96, fontweight='bold') 


y_min_ax1 = diversity_values_all.min()
y_max_ax1 = diversity_values_all.max()

# Set y-ticks to only the minimum and maximum values
ax1.set_yticks([y_min_ax1, y_max_ax1])

# Set y-tick labels with desired formatting
ax1.set_yticklabels([f"{y_min_ax1:.2f}", f"{y_max_ax1:.2f}"], fontsize=90, fontweight='bold')
# Align the x-axis of the line plot with the bar plot
ax1.set_xticks(range(len(sorted_loci_all)))
ax1.set_xticklabels(sorted_loci_all, rotation=45, ha='right', fontsize=96)
ax1.set_xlim([-0.5, len(sorted_loci_all) - 0.5])

# Ensure there is no x-axis label overlap with the lower plot
ax1.tick_params(axis='x', which='both', length=0)

# --- Bottom plot: Stacked bar plot ---
bottom_all = np.zeros(len(pivot_df_all.columns))  # Initialize the bottom of the bars

for allele in pivot_df_all.index:
    frequencies = pivot_df_all.loc[allele]
    # **Step 3: Map colors based on normalized global frequencies**
    colors = cmap(norm(frequencies.values))
    
    bars = ax2.bar(
        pivot_df_all.columns, 
        frequencies, 
        bottom=bottom_all, 
        color=colors,
        edgecolor='white',  # Add edge color for better visual separation
        label=allele
    )
    bottom_all += frequencies  # Update the bottom for the next allele

# Add title and labels
ax2.set_title('Allele Frequencies', fontsize=96, fontweight='bold')
# ax2.set_xlabel('Locus', fontsize=96)
ax2.set_ylabel('Allele Frequency', fontsize=96, fontweight='bold')


# Align the x-axis labels with the line plot
ax2.set_xticks(range(len(pivot_df_all.columns)))
ax2.set_xticklabels(pivot_df_all.columns, rotation=45, ha='right')
ax2.set_xlim([-0.5, len(pivot_df_all.columns) - 0.5])

ax2.tick_params(axis='y', which='both', labelsize=96)
ax2.tick_params(axis='x', which='both', labelsize=96)

# **Step 4: Create a color bar based on global frequency normalization**
sm_all = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm_all.set_array([])  # Set array to empty, as the colorbar is based on the normalized data

# Create the colorbar
cbar_all = plt.colorbar(sm_all, ax=ax2)

# Set the colorbar label with desired fontsize
cbar_all.set_label('Allele Frequency', fontsize=96, fontweight='bold')

# Set the tick labels fontsize
cbar_all.ax.tick_params(labelsize=96)  # Replace 40 with your desired fontsize

# Adjust layout for better spacing between plots
plt.tight_layout()

# Save the figure to an SVG file
plt.savefig(f'{outdir}/All_samples_KIR_allele_frequencies.pdf', format='pdf', bbox_inches='tight', dpi=50)

# Close the figure
plt.close()

print("Plot saved for all samples combined.")

# ---------------------------- Completion Message ----------------------------

print("All plots have been successfully generated and saved to the output directory.")

# Save the allele frequencies data to a CSV file for all samples
allele_counts_all[['Locus', 'Alleles', 'Frequency']].to_csv(f'{outdir}/allele_frequencies_all_samples.csv', index=False)
print(f"Allele frequencies saved to '{outdir}/allele_frequencies_all_samples.csv'.")