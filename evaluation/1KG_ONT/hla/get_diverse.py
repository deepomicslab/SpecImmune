import pandas as pd
import numpy as np
from tqdm import tqdm
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def get_population(df, sample_name):
    """
    Given a DataFrame and a sample name, return the population.
    
    Parameters:
    df (pd.DataFrame): The DataFrame containing the sample information.
    sample_name (str): The name of the sample to look up.
    
    Returns:
    str: The population of the sample.
    """
    # Filter the DataFrame to find the row with the given sample name
    sample_row = df[df['Sample'] == sample_name]
    
    # Check if the sample exists in the DataFrame
    if sample_row.empty:
        return "Unknown"
    
    # Extract and return the population
    return sample_row['Population'].values[0]

# Step 2: Read the gene read counts file into a DataFrame
gene_reads_df = pd.read_csv('../read_depth.csv', index_col=0)
# gene_reads_df.columns = gene_reads_df.columns.str.replace('HLA-', '', regex=False)

print(gene_reads_df)

# Step 3: Set the reads cutoff values
read_cutoffs = [15]  # Full list of cutoffs

# Step 5: Read the population and region data from the Excel file
metadata_df = pd.read_excel('./20130606_sample_info.xlsx', engine='openpyxl')

# Ensure consistent column names
metadata_df.columns = metadata_df.columns.str.replace(' ', '_')

allele_df = pd.read_csv('merged_samples.csv')

# Add population to allele_df
allele_df['Population'] = allele_df['Sample'].apply(lambda x: get_population(metadata_df, x))

# Split alleles in the Genotype column by ';' and only take the first allele
allele_df['Alleles'] = allele_df['Genotype'].str.split(';').str[0]

# Save the allele dataframe to a CSV file
allele_df.to_csv('allele_df.csv', index=False)

# Group by Population, Locus (gene), and Allele to count occurrences
allele_counts = allele_df.groupby(['Population', 'Locus', 'Alleles']).size().reset_index(name='Count')

# Add a "Total" column which is the sum of counts for the same Population and Locus
allele_counts['Total'] = allele_counts.groupby(['Population', 'Locus'])['Count'].transform('sum')

# Calculate the frequency of each allele by dividing the count by the total
allele_counts['Frequency'] = allele_counts['Count'] / allele_counts['Total']

# Save the allele counts dataframe to a CSV file
allele_counts.to_csv('allele_counts.csv', index=False)

# Calculate the Shannon diversity index for each Locus in each Population
def calculate_shannon_diversity(frequencies):
    """
    Calculate the Shannon diversity index for a set of frequencies.
    
    Parameters:
    frequencies (pd.Series): A series of allele frequencies.
    
    Returns:
    float: Shannon diversity index.
    """
    # Filter out zero frequencies to avoid log(0)
    frequencies = frequencies[frequencies > 0]
    return -np.sum(frequencies * np.log(frequencies))

# Group by Population and Locus, then calculate the diversity
diversity_df = allele_counts.groupby(['Population', 'Locus'])['Frequency'].apply(calculate_shannon_diversity).reset_index(name='Diversity')

# Save the diversity dataframe to a CSV file
diversity_df.to_csv('diversity_df.csv', index=False)

print(diversity_df)


# Assuming all previous steps (like reading data, calculating frequencies) are done
# and the final DataFrame `allele_counts` and `diversity_df` are ready

# Set the style for the plots
sns.set(style="white")

# # Get the list of unique populations
# populations = allele_counts['Population'].unique()

# # Loop through each population and generate the plot
# for population in populations:
#     # Filter data for the current population
#     pop_data = allele_counts[allele_counts['Population'] == population]

#     # Merge population data with diversity data
#     pop_data = pd.merge(pop_data, diversity_df, on=['Population', 'Locus'])

#     # Sort the data by Shannon diversity in descending order
#     pop_data = pop_data.sort_values(by='Diversity', ascending=False)

#     # Pivot the data to prepare for stacked bar plotting
#     pivot_df = pop_data.pivot(index='Alleles', columns='Locus', values='Frequency').fillna(0)

#     # Ensure the Locus columns are ordered by diversity
#     pivot_df = pivot_df[pop_data['Locus'].unique()]

#     # Sort alleles by their overall frequency for better visualization
#     pivot_df = pivot_df.reindex(pivot_df.sum(axis=1).sort_values(ascending=False).index)

#     # Create a figure with two subplots: one for the line plot and one for the stacked bar plot
#     fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(10, 5), gridspec_kw={'height_ratios': [1, 3]})

#     # --- Top plot: Shannon diversity line plot ---
#     # Adjust the aspect ratio to control the 'apparent' width
#     ax1.set_aspect(aspect='auto')

#     # Plot Shannon diversity for each Locus, using the pivot_df column order
#     diversity_values = pop_data.groupby('Locus')['Diversity'].first().reindex(pivot_df.columns)
#     ax1.plot(pivot_df.columns, diversity_values, marker='o', color='red')

#     # Add title and labels
#     ax1.set_title('Shannon Diversity Across Loci', fontsize=16)
#     ax1.set_ylabel('Diversity', fontsize=14)

#     # Align the x-axis of the line plot with the bar plot
#     ax1.set_xticks(range(len(pivot_df.columns)))
#     ax1.set_xticklabels(pivot_df.columns, rotation=45, ha='right')
#     ax1.set_xlim([-0.5, len(pivot_df.columns) - 0.5])

#     # Ensure there is no x-axis label overlap with the lower plot
#     ax1.tick_params(axis='x', which='both', length=0)

#     # --- Bottom plot: Stacked bar plot ---

#     # Define a colormap
#     cmap = plt.cm.viridis

#     # Normalize the frequency values for color mapping
#     norm = plt.Normalize(0, 1)

#     # Plot each allele with its corresponding color
#     bottom = np.zeros(len(pivot_df.columns))  # Initialize the bottom of the bars
#     for allele in pivot_df.index:
#         frequencies = pivot_df.loc[allele]
#         bars = ax2.bar(
#             pivot_df.columns, 
#             frequencies, 
#             bottom=bottom, 
#             color=cmap(norm(frequencies.values)),
#             edgecolor='white',  # Add edge color for better visual separation
#             label=allele
#         )
#         bottom += frequencies  # Update the bottom for the next allele

#     # Add title and labels
#     ax2.set_title(f'Population: {population}', fontsize=16)
#     ax2.set_xlabel('Gene (Locus)', fontsize=14)
#     ax2.set_ylabel('Allele Frequency', fontsize=14)

#     # Align the x-axis labels with the line plot
#     ax2.set_xticks(range(len(pivot_df.columns)))
#     ax2.set_xticklabels(pivot_df.columns, rotation=45, ha='right')
#     ax2.set_xlim([-0.5, len(pivot_df.columns) - 0.5])

#     # Create a color bar
#     sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
#     sm.set_array([])  # Set array to empty, as the color bar is based on the normalized data
#     cbar = plt.colorbar(sm, ax=ax2)
#     cbar.set_label('Frequency', fontsize=14)

#     # Adjust layout for better spacing between plots
#     plt.tight_layout()

#     # Save the figure to a PDF file
#     # plt.savefig(f'{population}_allele_frequencies.pdf')
#     # plt.savefig(f'{population}_allele_frequencies.pdf', format='pdf', bbox_inches='tight')
#     plt.savefig(f'{population}_allele_frequencies.png', format='png', bbox_inches='tight')
#     # Optionally close the figure to free up memory if looping through many populations
#     plt.close()

# After reading the allele_df and before counting populations, add:

# Remove duplicates based on the 'Sample' column
allele_df_unique = allele_df.drop_duplicates(subset='Sample')

# Now use this deduplicated DataFrame for population counting
population_counts = allele_df_unique['Population'].value_counts().reset_index()
population_counts.columns = ['Population', 'Count']

# Sort the populations by count in descending order
population_counts = population_counts.sort_values('Count', ascending=False)

# Save the population counts to a CSV file
population_counts.to_csv('population_counts.csv', index=False)

print("Population counts (after removing duplicates):")
print(population_counts)

# Create a bar plot of the population counts
plt.figure(figsize=(12, 6))
sns.barplot(x='Population', y='Count', data=population_counts)
plt.title('Number of Unique Samples per Population')
plt.xlabel('Population')
plt.ylabel('Number of Unique Samples')
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.savefig('population_counts.png')
plt.close()

# Print the total number of unique samples
total_unique_samples = population_counts['Count'].sum()
print(f"\nTotal number of unique samples: {total_unique_samples}")