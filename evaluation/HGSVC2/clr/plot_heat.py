import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import math

# Step 1: Read the CSV files into DataFrames
df1 = pd.read_csv('hlala.match.csv')
df2 = pd.read_csv('spechla.match.csv')
df3 = pd.read_csv('speclong.match.csv')

# Combine the DataFrames into a list
dfs = [df1, df2, df3]
df_names = ['HLA*LA', 'SpecHLA', 'SpecLong']  # Updated labels for DataFrames

# Define colors to be used for each software tool
colors = ["#FF5959", "#FFAD5A", "#4F9DA6"]  # Colors for HLA*LA, SpecHLA, and SpecLong respectively

# Step 2: Read the gene read counts file into a DataFrame
gene_reads_df = pd.read_csv('read_depth.csv', index_col=0)
gene_reads_df.columns = gene_reads_df.columns.str.replace('HLA-', '', regex=False)

# Define the specific genes of interest
target_genes = ['A', 'B', 'C', 'DQB1', 'DQA1']

# Filter gene_reads_df to only include the target genes
gene_reads_df = gene_reads_df[target_genes]

# Step 3: Set the reads cutoff values
read_cutoffs = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50]  # Full list of cutoffs

# Step 4: Initialize a dictionary to store accuracy results for plotting
plot_data = {cutoff: {df_name: [] for df_name in df_names} for cutoff in read_cutoffs}

# Calculate detailed accuracy for each gene in each DataFrame based on each reads cutoff
for cutoff in read_cutoffs:
    # Filter gene_reads_df based on the cutoff (samples below cutoff will be NaN)
    filtered_gene_reads_df = gene_reads_df.where(gene_reads_df >= cutoff)

    # Loop through each DataFrame and calculate accuracy
    for df, df_name in zip(dfs, df_names):
        # Calculate accuracy for each gene
        for gene in target_genes:
            # Get the samples that pass the filter for this gene
            valid_samples = filtered_gene_reads_df[gene].dropna().index

            # Filter the current DataFrame to consider only the valid samples for this gene
            filtered_df = df[(df['Gene'] == gene) & (df['Sample'].isin(valid_samples))]

            # If there are no valid samples, skip to the next gene
            if filtered_df.empty:
                continue

            # Sum the Match and Total for the valid samples
            total_match = filtered_df['Match'].sum()
            total_total = filtered_df['Total'].sum()

            # Calculate accuracy for the gene
            accuracy = total_match / total_total if total_total > 0 else None

            # Store the result in the plot_data dictionary
            plot_data[cutoff][df_name].append({'Gene': gene, 'Accuracy': accuracy})

# Now, plot_data contains the accuracy for the specified genes
print(plot_data)