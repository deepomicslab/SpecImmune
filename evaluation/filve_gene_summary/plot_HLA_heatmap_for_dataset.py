import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# Define your datasets with corresponding names
datasets = ['../1KG_ONT/hla/', '../HPRC/hla/hifi/', '../HPRC/hla/ont/', '../HGSVC2/hifi', '../HGSVC2/clr']
dataset_names = ['1KG ONT', 'HPRC HiFi', 'HPRC ONT', 'HGSVC2 HiFi', 'HGSVC2 CLR']
csv_files = ['hlala.match.csv', 'spechla.match.csv', 'speclong.match.csv']
df_names = ['HLA*LA', 'SpecHLA', 'SpecLong']
target_genes = ['A', 'B', 'C', 'DQB1', 'DRB1']

# Set read cutoff
read_cutoff = 20

# Initialize a dictionary to store the accuracy results
accuracy_data = []

for dataset, dataset_name in zip(datasets, dataset_names):
    # Read the gene read counts for the current dataset
    gene_reads_df = pd.read_csv(f'{dataset}/read_depth.csv', index_col=0)
    gene_reads_df.columns = gene_reads_df.columns.str.replace('HLA-', '', regex=False)

    # Filter gene_reads_df to only include the target genes
    gene_reads_df = gene_reads_df[target_genes]

    # Apply the read cutoff
    filtered_gene_reads_df = gene_reads_df.where(gene_reads_df >= read_cutoff)

    for csv_file, df_name in zip(csv_files, df_names):
        # Read the match data for the current software tool and dataset
        df = pd.read_csv(f'{dataset}/{csv_file}')

        for gene in target_genes:
            # Get the samples that pass the filter for this gene
            valid_samples = filtered_gene_reads_df[gene].dropna().index

            # Filter the current DataFrame to consider only the valid samples for this gene
            filtered_df = df[(df['Gene'] == gene) & (df['Sample'].isin(valid_samples))]

            # Calculate accuracy for the gene
            accuracy = None
            if not filtered_df.empty:
                total_match = filtered_df['Match'].sum()
                total_total = filtered_df['Total'].sum()
                accuracy = total_match / total_total if total_total > 0 else None

            # Store the result in the accuracy_data list
            accuracy_data.append({
                'Dataset': dataset_name,  # Use the simplified dataset name
                'Software': df_name,
                'Gene': gene,
                'Accuracy': accuracy
            })

# Convert the accuracy data to a DataFrame
accuracy_df = pd.DataFrame(accuracy_data)

# Define custom colors
colors = ["#FF5959", "#FFAD5A", "#4F9DA6"]
colors = ["#D2D2D2", "#062565", "#0098B4"]

# Plot and save grouped bar plots for each dataset
for dataset_name in accuracy_df['Dataset'].unique():
    plt.figure(figsize=(6, 3))
    
    # Filter data for the current dataset
    dataset_df = accuracy_df[accuracy_df['Dataset'] == dataset_name]
    
    # Create a barplot with custom colors
    sns.barplot(
        x='Gene', 
        y='Accuracy', 
        hue='Software', 
        data=dataset_df, 
        palette=colors, 
        ci=None
    )
    
    plt.title(f'{dataset_name} - Gene Accuracy', fontsize=16)
    plt.ylim(0, 1.1)
    plt.xticks(rotation=45)
    plt.ylabel('Accuracy', fontsize=14)
    plt.xlabel('Locus', fontsize=14)
    # set label font size
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    
    # Adjusting legend placement
    plt.legend( 
        loc='center left', 
        bbox_to_anchor=(1, 0.5), 
        fontsize=14, 
        frameon=False
    )
    
    plt.tight_layout()
    
    # Save the plot to a PDF file using the simplified dataset name
    plt.savefig(f"{dataset_name.replace(' ', '_')}_5_gene_accuracy.svg", format='svg', bbox_inches='tight', dpi=600)
    
    # Close the plot to free up memory
    plt.close()