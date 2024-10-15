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
colors = ["#D2D2D2", "#062565", "#0098B4"]


# Step 2: Read the gene read counts file into a DataFrame
gene_reads_df = pd.read_csv('read_depth.csv', index_col=0)
gene_reads_df.columns = gene_reads_df.columns.str.replace('HLA-', '', regex=False)

# Define gene classes
class1_genes = ['A', 'B', 'C', 'E', 'F', 'G']
# class1_pseudogenes = ['H', 'J', 'K', 'L', 'N', 'P', 'R', 'S', 'T', 'U', 'V', 'W', 'Y']
class1_pseudogenes = ['H', 'J', 'K', 'L', 'N', 'P', 'S', 'T', 'U', 'V', 'W']

# class2_genes = ['DRA', 'DRB', 'DQA1', 'DQA2', 'DQB1', 'DQB2', 'DPA1', 'DPA2', 'DPB1', 'DPB2', 'DMA', 'DMB', 'DOA', 'DOB', 
#                 'DRB1', 'DRB2', 'DRB3', 'DRB4', 'DRB5', 'DRB6', 'DRB7', 'DRB8', 'DRB9']

class2_genes = ['DRA', 'DQA1', 'DQA2', 'DQB1', 'DPA1', 'DPA2', 'DPB1', 'DPB2', 'DMA', 'DMB', 'DOA', 'DOB', 
                'DRB1', 'DRB3', 'DRB4']
# non_hla_genes = ['MICA', 'MICB', 'TAP1', 'TAP2', 'HFE']
non_hla_genes = ['MICA', 'MICB', 'TAP1', 'TAP2']


# Filter gene_reads_df to only include genes present in any of the DataFrames
genes_in_dfs = pd.concat([df['Gene'] for df in dfs]).unique()
gene_reads_df = gene_reads_df[genes_in_dfs]

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
        for gene in genes_in_dfs:
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

# Step 5: Plot the performance of the software across different gene classes and cutoffs

# Figure 1: Gene classes vs. cutoff
plt.figure(figsize=(7.0, 5.0))  # Adjusted figure size for double-column width

# Define gene classes and their corresponding genes
gene_classes = {
    'Class I': class1_genes,
    'Class II': class2_genes,
    'Pseudogenes': class1_pseudogenes,
    'Non-HLA': non_hla_genes
}
class_colors = {
    'Class I': '#1f77b4',
    'Class II': '#ff7f0e',
    'Pseudogenes': '#2ca02c',
    'Non-HLA': '#d62728'
}

# Create a 2x2 grid of subplots for the gene classes
for i, (class_name, genes) in enumerate(gene_classes.items(), 1):
    plt.subplot(2, 2, i)  # 2x2 grid of subplots

    for j, (df_name, color) in enumerate(zip(df_names, colors)):
        accuracies_by_cutoff = []
        for cutoff in read_cutoffs:
            # Calculate the average accuracy for the genes in the current class
            gene_accuracies = [item['Accuracy'] for item in plot_data[cutoff][df_name] if item['Gene'] in genes]
            avg_accuracy = sum(gene_accuracies) / len(gene_accuracies) if gene_accuracies else None
            accuracies_by_cutoff.append(avg_accuracy)

        sns.lineplot(x=read_cutoffs, y=accuracies_by_cutoff, label=df_name, color=color, marker='o')

    plt.title(f"{class_name}", fontsize=14)
    plt.xlabel("Depth Cutoff", fontsize=14)
    plt.ylabel("Average Accuracy", fontsize=10)
    plt.xticks(rotation=45, fontsize=14)
    plt.yticks(fontsize=14)
    plt.ylim(0, 1.1)  # Assuming accuracy is between 0 and 1

    # Legend inside the plot
    plt.legend(title="", fontsize=8, frameon=False)

plt.tight_layout(pad=1)  # Adjust padding for better spacing
plt.savefig("HPRC_HIFI_performance_by_class_and_cutoff.pdf", bbox_inches='tight', format="pdf", dpi=600)  # Save the plot as a PDF

# Figure 2: Individual gene plots grouped by class, with 4 subplots per row
for class_name, genes in gene_classes.items():
    rows = math.ceil(len(genes) / 4)  # Calculate the number of rows needed
    
    plt.figure(figsize=(14.0, rows * 3))  # Consistent figure size: 14 inches wide and 3 inches per row

    for idx, gene in enumerate(genes):
        plt.subplot(rows, 4, idx + 1)  # 4 subplots per row

        for df_name, color in zip(df_names, colors):
            accuracies_by_cutoff = []
            for cutoff in read_cutoffs:
                # Retrieve accuracy data for the gene
                gene_accuracies = [item['Accuracy'] for item in plot_data[cutoff][df_name] if item['Gene'] == gene]
                avg_accuracy = gene_accuracies[0] if gene_accuracies else None
                accuracies_by_cutoff.append(avg_accuracy)
            
            sns.lineplot(x=read_cutoffs, y=accuracies_by_cutoff, label=df_name, color=color, marker='o')

        plt.title(f"{gene}", fontsize=14)
        plt.xlabel("Depth Cutoff", fontsize=14)
        plt.ylabel("Accuracy", fontsize=14)
        plt.xticks(rotation=45, fontsize=14)
        plt.yticks(fontsize=14)
        plt.ylim(0, 1.1)

        # Add the legend only to the first subplot in the first row
        # if idx == 0:
        plt.legend(title="Methods", fontsize=8, frameon=False, loc="lower left")

    plt.tight_layout(pad=1)  # Adjust padding for better spacing
    plt.savefig(f"HPRC_HIFI_performance_by_gene_and_cutoff_{class_name}.pdf", bbox_inches='tight', format="pdf", dpi=600)  # Save the plot as a PDF

# Uncomment the line below if you want to display the plot interactively
# plt.show()