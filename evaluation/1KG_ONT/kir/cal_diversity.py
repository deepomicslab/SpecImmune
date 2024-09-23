import pandas as pd
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
import seaborn as sns

# Create a population to superpopulation mapping dictionary
population_superpopulation_dict = {
    'CDX': 'EAS',  # Chinese Dai in Xishuangbanna, China
    'CHB': 'EAS',  # Han Chinese in Bejing, China
    'JPT': 'EAS',  # Japanese in Tokyo, Japan
    'KHV': 'EAS',  # Kinh in Ho Chi Minh City, Vietnam
    'CHS': 'EAS',  # Southern Han Chinese, China
    'BEB': 'SAS',  # Bengali in Bangladesh
    'GIH': 'SAS',  # Gujarati Indian in Houston,TX
    'ITU': 'SAS',  # Indian Telugu in the UK
    'PJL': 'SAS',  # Punjabi in Lahore,Pakistan
    'STU': 'SAS',  # Sri Lankan Tamil in the UK
    'ASW': 'AFR',  # African Ancestry in Southwest US
    'ACB': 'AFR',  # African Caribbean in Barbados
    'ESN': 'AFR',  # Esan in Nigeria
    'GWD': 'AFR',  # Gambian in Western Division, The Gambia
    'LWK': 'AFR',  # Luhya in Webuye, Kenya
    'MSL': 'AFR',  # Mende in Sierra Leone
    'YRI': 'AFR',  # Yoruba in Ibadan, Nigeria
    'GBR': 'EUR',  # British in England and Scotland
    'FIN': 'EUR',  # Finnish in Finland
    'IBS': 'EUR',  # Iberian populations in Spain
    'TSI': 'EUR',  # Toscani in Italy
    'CEU': 'EUR',  # Utah residents with Northern and Western European ancestry
    'CLM': 'AMR',  # Colombian in Medellin, Colombia
    'MXL': 'AMR',  # Mexican Ancestry in Los Angeles, California
    'PEL': 'AMR',  # Peruvian in Lima, Peru
    'PUR': 'AMR'   # Puerto Rican in Puerto Rico
}

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

# Step 2: Read the gene read counts file into a DataFrame
gene_reads_df = pd.read_csv('./read_depth_kirs.csv', index_col=0)

# Step 3: Set the reads cutoff values
read_cutoffs = [10]  # Full list of cutoffs

# Step 5: Read the population and region data from the Excel file
metadata_df = pd.read_excel('./20130606_sample_info.xlsx', engine='openpyxl')

# Ensure consistent column names
metadata_df.columns = metadata_df.columns.str.replace(' ', '_')

# Step 6: Read allele genotype data
allele_df = pd.read_csv('merged_samples.csv')

# Add population to allele_df
allele_df['Population'] = allele_df['Sample'].apply(lambda x: get_population(metadata_df, x))

# Split alleles in the Genotype column by ';' and only take the first allele
allele_df['Alleles'] = allele_df['Genotype'].str.split(';').str[0]

# Save the allele dataframe to a CSV file
allele_df.to_csv('allele_df.csv', index=False)

# Step 7: Group by Population, Locus (gene), and Alleles to count occurrences
allele_counts = allele_df.groupby(['Population', 'Locus', 'Alleles']).size().reset_index(name='Count')

# Add a "Total" column which is the sum of counts for the same Population and Locus
allele_counts['Total'] = allele_counts.groupby(['Population', 'Locus'])['Count'].transform('sum')

# Calculate the frequency of each allele by dividing the count by the total
allele_counts['Frequency'] = allele_counts['Count'] / allele_counts['Total']

# Save the allele counts dataframe to a CSV file
allele_counts.to_csv('allele_counts.csv', index=False)

# Step 8: Define Shannon diversity function
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

# Step 9: Group by Population and Locus, then calculate the diversity
diversity_df = allele_counts.groupby(['Population', 'Locus'])['Frequency'].apply(calculate_shannon_diversity).reset_index(name='Diversity')

# Save the diversity dataframe to a CSV file
diversity_df.to_csv('diversity_df.csv', index=False)

print("Population-level diversity:")
print(diversity_df)

# Step 10: Calculate diversity at the superpopulation level

# Map each population to its superpopulation
allele_counts['Superpopulation'] = allele_counts['Population'].map(population_superpopulation_dict)

# Group by Superpopulation, Locus, and Alleles to aggregate counts
superpopulation_allele_counts = allele_counts.groupby(['Superpopulation', 'Locus', 'Alleles'])['Count'].sum().reset_index()

# Calculate the total count for each Superpopulation and Locus
superpopulation_allele_counts['Total'] = superpopulation_allele_counts.groupby(['Superpopulation', 'Locus'])['Count'].transform('sum')

# Calculate the frequency of each allele at the superpopulation level
superpopulation_allele_counts['Frequency'] = superpopulation_allele_counts['Count'] / superpopulation_allele_counts['Total']

# Group by Superpopulation and Locus, then calculate the diversity
superpopulation_diversity_df = superpopulation_allele_counts.groupby(['Superpopulation', 'Locus'])['Frequency'].apply(calculate_shannon_diversity).reset_index(name='Diversity')

# Save the superpopulation diversity dataframe to a CSV file
superpopulation_diversity_df.to_csv('superpopulation_diversity_df_KIR.csv', index=False)

print("Superpopulation-level diversity:")
print(superpopulation_diversity_df)


