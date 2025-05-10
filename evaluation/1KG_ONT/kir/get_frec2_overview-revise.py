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

def calculate_zahl_jackknife_exact(counts):
    """
    Calculate the Zahl-Jackknife exact method for entropy estimation.
    
    Parameters:
    counts (array-like): Array of counts.
    
    Returns:
    float: Estimated diversity index.
    """
    counts = np.array(counts)
    counts = counts[counts > 0]
    n = len(counts)
    total = counts.sum()
    freqs = counts / total
    H_full = -np.sum(freqs * np.log(freqs))
    
    if n <= 1:  # Can't jackknife with only one element
        return H_full
    
    H_jack_list = []
    for i in range(n):
        reduced_counts = np.delete(counts, i)
        reduced_total = reduced_counts.sum()
        reduced_freqs = reduced_counts / reduced_total
        H_i = -np.sum(reduced_freqs * np.log(reduced_freqs))
        H_jack_list.append(H_i)
    H_jack_mean = np.mean(H_jack_list)
    return n * H_full - (n - 1) * H_jack_mean

def calculate_chao_entropy_exact(counts):
    """
    Calculate the Chao entropy exact method for diversity estimation.
    
    Parameters:
    counts (array-like): Array of counts.
    
    Returns:
    float: Estimated diversity index.
    """
    counts = np.array(counts)
    counts = counts[counts > 0]
    n = counts.sum()
    freqs = counts / n
    H_obs = -np.sum(freqs * np.log(freqs))
    
    f1 = np.sum(counts == 1)
    f2 = np.sum(counts == 2)
    
    if f1 == 0:
        return H_obs
    
    if f2 == 0:
        correction = (f1 * (f1 - 1)) / (2 * (f2 + 1))
    else:
        correction = (f1 ** 2) / (2 * f2)
    
    A = correction
    bias_term = (f1 / n) * np.log(n / f1) + A / n
    return H_obs + bias_term

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
allele_df.to_csv(f'{outdir}/allele_df.csv', index=False)
print("Allele DataFrame saved to 'allele_df.csv'.")

# Group by Population, Locus (gene), and Allele to count occurrences
allele_counts = allele_df.groupby(['Population', 'Locus', 'Alleles']).size().reset_index(name='Count')

# Add a "Total" column which is the sum of counts for the same Population and Locus
allele_counts['Total'] = allele_counts.groupby(['Population', 'Locus'])['Count'].transform('sum')

# Calculate the frequency of each allele by dividing the count by the total
allele_counts['Frequency'] = allele_counts['Count'] / allele_counts['Total']

# Save the allele counts DataFrame to a CSV file
allele_counts.to_csv(f'{outdir}/allele_counts.csv', index=False)
print("Allele counts saved to 'allele_counts.csv'.")

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

# ---------------------------- Calculate All Three Diversity Measures ----------------------------

# Calculate Shannon diversity for all samples per locus
diversity_shannon = allele_counts_all.groupby('Locus')['Frequency'].apply(calculate_shannon_diversity).reset_index(name='Diversity')
print("Shannon diversity calculation completed.")

# Calculate Zahl-Jackknife diversity for all samples per locus
diversity_jackknife = pd.DataFrame()
for locus, group in allele_counts_all.groupby('Locus'):
    jackknife = calculate_zahl_jackknife_exact(group['Count'].values)
    diversity_jackknife = pd.concat([diversity_jackknife, pd.DataFrame([{'Locus': locus, 'Diversity': jackknife}])], ignore_index=True)
print("Zahl-Jackknife diversity calculation completed.")

# Calculate Chao diversity for all samples per locus
diversity_chao = pd.DataFrame()
for locus, group in allele_counts_all.groupby('Locus'):
    chao = calculate_chao_entropy_exact(group['Count'].values)
    diversity_chao = pd.concat([diversity_chao, pd.DataFrame([{'Locus': locus, 'Diversity': chao}])], ignore_index=True)
print("Chao diversity calculation completed.")

# Save all three diversity measures to CSV files
diversity_shannon.to_csv(f'{outdir}/diversity_shannon.csv', index=False)
diversity_jackknife.to_csv(f'{outdir}/diversity_jackknife.csv', index=False)
diversity_chao.to_csv(f'{outdir}/diversity_chao.csv', index=False)

# Sort each diversity dataframe by diversity value in descending order
diversity_shannon_sorted = diversity_shannon.sort_values(by='Diversity', ascending=False)
diversity_jackknife_sorted = diversity_jackknife.sort_values(by='Diversity', ascending=False)
diversity_chao_sorted = diversity_chao.sort_values(by='Diversity', ascending=False)

# Extract the sorted loci for each method
sorted_loci_shannon = diversity_shannon_sorted['Locus'].tolist()
sorted_loci_jackknife = diversity_jackknife_sorted['Locus'].tolist()
sorted_loci_chao = diversity_chao_sorted['Locus'].tolist()

# Pivot the data for plotting
pivot_df_all = allele_counts_all.pivot(index='Alleles', columns='Locus', values='Frequency').fillna(0)

# Sort alleles by their overall frequency for better visualization
pivot_df_all = pivot_df_all.reindex(pivot_df_all.sum(axis=1).sort_values(ascending=False).index)

# ---------------------------- Define Global Normalization Based on Frequencies ----------------------------

# Calculate global min and max frequencies across all alleles and loci
global_min_frequency = allele_counts_all['Frequency'].min()
global_max_frequency = allele_counts_all['Frequency'].max()

print(f"Global Minimum Allele Frequency: {global_min_frequency}")
print(f"Global Maximum Allele Frequency: {global_max_frequency}")

# Define normalization based on global min and max frequencies
norm = plt.Normalize(vmin=global_min_frequency, vmax=global_max_frequency)

# ---------------------------- Function to Create Plot for a Specific Diversity Method ----------------------------

def create_diversity_plot(diversity_df, sorted_loci, method_name, outdir):
    """
    Create a plot showing diversity values and allele frequencies for a specific diversity method.
    
    Parameters:
    diversity_df (pd.DataFrame): DataFrame containing diversity values for each locus
    sorted_loci (list): List of loci sorted by diversity in descending order
    method_name (str): Name of the diversity method being plotted
    outdir (str): Directory to save output files
    """
    # Reorder pivot dataframe columns based on the sorted loci
    pivot_df_method = pivot_df_all[sorted_loci].copy()
    
    # Plot diversity measure alongside allele frequencies
    fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(120, 90), gridspec_kw={'height_ratios': [1, 10]})
    
    # --- Top plot: Diversity line plot ---
    diversity_values = diversity_df.set_index('Locus').loc[sorted_loci, 'Diversity'].values
    ax1.plot(sorted_loci, diversity_values, marker='o', color='red', linewidth=6, markersize=30)
    
    # Add title and labels
    ax1.set_title(f'{method_name} Diversity across KIR genes', fontsize=96, fontweight='bold')
    ax1.set_ylabel('Diversity', fontsize=96, fontweight='bold')
    
    # Set y-ticks to only the minimum and maximum values
    y_min_ax1 = diversity_values.min()
    y_max_ax1 = diversity_values.max()
    ax1.set_yticks([y_min_ax1, y_max_ax1])
    ax1.set_yticklabels([f"{y_min_ax1:.2f}", f"{y_max_ax1:.2f}"], fontsize=90, fontweight='bold')
    
    # Align the x-axis of the line plot with the bar plot
    ax1.set_xticks(range(len(sorted_loci)))
    ax1.set_xticklabels(sorted_loci, rotation=45, ha='right', fontsize=96)
    ax1.set_xlim([-0.5, len(sorted_loci) - 0.5])
    
    # Ensure there is no x-axis label overlap with the lower plot
    ax1.tick_params(axis='x', which='both', length=0)
    
    # --- Bottom plot: Stacked bar plot ---
    bottom_all = np.zeros(len(pivot_df_method.columns))  # Initialize the bottom of the bars
    
    for allele in pivot_df_method.index:
        frequencies = pivot_df_method.loc[allele]
        # Map colors based on normalized global frequencies
        colors = cmap(norm(frequencies.values))
        
        bars = ax2.bar(
            pivot_df_method.columns, 
            frequencies, 
            bottom=bottom_all, 
            color=colors,
            edgecolor='white',  # Add edge color for better visual separation
            label=allele
        )
        bottom_all += frequencies  # Update the bottom for the next allele
    
    # Add title and labels
    ax2.set_title('Allele Frequencies', fontsize=96, fontweight='bold')
    ax2.set_ylabel('Allele Frequency', fontsize=96, fontweight='bold')
    
    # Align the x-axis labels with the line plot
    ax2.set_xticks(range(len(pivot_df_method.columns)))
    ax2.set_xticklabels(pivot_df_method.columns, rotation=45, ha='right')
    ax2.set_xlim([-0.5, len(pivot_df_method.columns) - 0.5])
    
    ax2.tick_params(axis='y', which='both', labelsize=96)
    ax2.tick_params(axis='x', which='both', labelsize=96)
    
    # Create a color bar based on global frequency normalization
    sm_all = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm_all.set_array([])  # Set array to empty, as the colorbar is based on the normalized data
    
    # Create the colorbar
    cbar_all = plt.colorbar(sm_all, ax=ax2)
    
    # Set the colorbar label with desired fontsize
    cbar_all.set_label('Allele Frequency', fontsize=96, fontweight='bold')
    
    # Set the tick labels fontsize
    cbar_all.ax.tick_params(labelsize=96)
    
    # Adjust layout for better spacing between plots
    plt.tight_layout()
    
    # Save the figure to PDF and PNG files
    plt.savefig(f'{outdir}/All_samples_KIR_allele_frequencies_{method_name}.pdf', format='pdf', bbox_inches='tight', dpi=50)
    plt.savefig(f'{outdir}/All_samples_KIR_allele_frequencies_{method_name}.png', bbox_inches='tight', dpi=150)
    
    # Close the figure
    plt.close()
    
    print(f"Plot saved for all samples combined with {method_name} diversity.")

# ---------------------------- Create Plots for Each Diversity Method ----------------------------

# Create plots for each diversity method
create_diversity_plot(diversity_shannon_sorted, sorted_loci_shannon, 'Shannon', outdir)
create_diversity_plot(diversity_jackknife_sorted, sorted_loci_jackknife, 'Jackknife', outdir)
create_diversity_plot(diversity_chao_sorted, sorted_loci_chao, 'Chao', outdir)

# ---------------------------- Original Shannon Plot (For Compatibility) ----------------------------

# This creates the original format plot with Shannon diversity for compatibility
fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(120, 90), gridspec_kw={'height_ratios': [1, 10]})

# Reorder pivot dataframe columns based on Shannon diversity
pivot_df_shannon = pivot_df_all[sorted_loci_shannon]

# --- Top plot: Shannon diversity line plot ---
diversity_values_shannon = diversity_shannon_sorted['Diversity'].values
ax1.plot(sorted_loci_shannon, diversity_values_shannon, marker='o', color='red', linewidth=6, markersize=30)

# Add title and labels
ax1.set_title('Shannon diversity across KIR genes', fontsize=96, fontweight='bold')
ax1.set_ylabel('Diversity', fontsize=96, fontweight='bold')

y_min_ax1 = diversity_values_shannon.min()
y_max_ax1 = diversity_values_shannon.max()

# Set y-ticks to only the minimum and maximum values
ax1.set_yticks([y_min_ax1, y_max_ax1])

# Set y-tick labels with desired formatting
ax1.set_yticklabels([f"{y_min_ax1:.2f}", f"{y_max_ax1:.2f}"], fontsize=90, fontweight='bold')
# Align the x-axis of the line plot with the bar plot
ax1.set_xticks(range(len(sorted_loci_shannon)))
ax1.set_xticklabels(sorted_loci_shannon, rotation=45, ha='right', fontsize=96)
ax1.set_xlim([-0.5, len(sorted_loci_shannon) - 0.5])

# Ensure there is no x-axis label overlap with the lower plot
ax1.tick_params(axis='x', which='both', length=0)

# --- Bottom plot: Stacked bar plot ---
bottom_all = np.zeros(len(pivot_df_shannon.columns))  # Initialize the bottom of the bars

for allele in pivot_df_shannon.index:
    frequencies = pivot_df_shannon.loc[allele]
    # Map colors based on normalized global frequencies
    colors = cmap(norm(frequencies.values))
    
    bars = ax2.bar(
        pivot_df_shannon.columns, 
        frequencies, 
        bottom=bottom_all, 
        color=colors,
        edgecolor='white',  # Add edge color for better visual separation
        label=allele
    )
    bottom_all += frequencies  # Update the bottom for the next allele

# Add title and labels
ax2.set_title('Allele Frequencies', fontsize=96, fontweight='bold')
ax2.set_ylabel('Allele Frequency', fontsize=96, fontweight='bold')

# Align the x-axis labels with the line plot
ax2.set_xticks(range(len(pivot_df_shannon.columns)))
ax2.set_xticklabels(pivot_df_shannon.columns, rotation=45, ha='right')
ax2.set_xlim([-0.5, len(pivot_df_shannon.columns) - 0.5])

ax2.tick_params(axis='y', which='both', labelsize=96)
ax2.tick_params(axis='x', which='both', labelsize=96)

# Create a color bar based on global frequency normalization
sm_all = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm_all.set_array([])  # Set array to empty, as the colorbar is based on the normalized data

# Create the colorbar
cbar_all = plt.colorbar(sm_all, ax=ax2)

# Set the colorbar label with desired fontsize
cbar_all.set_label('Allele Frequency', fontsize=96, fontweight='bold')

# Set the tick labels fontsize
cbar_all.ax.tick_params(labelsize=96)

# Adjust layout for better spacing between plots
plt.tight_layout()

# Save the figure
plt.savefig(f'{outdir}/All_samples_KIR_allele_frequencies.pdf', format='pdf', bbox_inches='tight', dpi=50)

# Close the figure
plt.close()

print("Original format plot with Shannon diversity saved.")

# ---------------------------- Save Allele Frequencies Data for All Samples ----------------------------

# Save the allele frequencies data to a CSV file for all samples
allele_counts_all[['Locus', 'Alleles', 'Frequency']].to_csv(f'{outdir}/allele_frequencies_all_samples.csv', index=False)
print(f"Allele frequencies saved to '{outdir}/allele_frequencies_all_samples.csv'.")

# ---------------------------- Completion Message ----------------------------

print("All plots have been successfully generated and saved to the output directory.")