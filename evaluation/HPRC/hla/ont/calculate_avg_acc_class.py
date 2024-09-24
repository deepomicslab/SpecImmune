import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Step 1: Read the CSV files into DataFrames
df1 = pd.read_csv('hlala.match.csv')
df2 = pd.read_csv('spechla.match.csv')
df3 = pd.read_csv('speclong.match.csv')

# Combine the DataFrames into a list
dfs = [df1, df2, df3]
df_names = ['HLA*LA', 'SpecHLA', 'SpecLong']  # Updated labels for DataFrames

# Step 2: Define the genes to consider for each software
genes_for_hlala = ['A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRA', 'DRB1', 'E', 'F', 'G', 'H', 'K', 'V', 'DRB3', 'DRB4']
genes_for_spechla = ['A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRB1']
genes_for_speclong = [gene for gene in pd.concat([df1['Gene'], df2['Gene'], df3['Gene']]).unique() if gene not in ['HFE']]

# Step 3: Filter the gene read counts to only include relevant genes for each software
gene_reads_df = pd.read_csv('read_depth.csv', index_col=0)
gene_reads_df.columns = gene_reads_df.columns.str.replace('HLA-', '', regex=False)

# Define a dictionary to map software names to their respective gene lists
software_genes = {
    'HLA*LA': genes_for_hlala,
    'SpecHLA': genes_for_spechla,
    'SpecLong': genes_for_speclong
}

# Define gene classes
class1_genes = ['A', 'B', 'C', 'E', 'F', 'G']
class1_pseudogenes = ['H', 'J', 'K', 'L', 'N', 'P', 'S', 'T', 'U', 'V', 'W']
class2_genes = ['DRA', 'DQA1', 'DQA2', 'DQB1', 'DPA1', 'DPA2', 'DPB1', 'DPB2', 'DMA', 'DMB', 'DOA', 'DOB', 'DRB1', 'DRB3', 'DRB4']
non_hla_genes = ['MICA', 'MICB', 'TAP1', 'TAP2']

gene_classes = {
    'Class I Genes': class1_genes,
    'Class I Pseudogenes': class1_pseudogenes,
    'Class II Genes': class2_genes,
    'Non-HLA Genes': non_hla_genes
}

# Step 4: Set the chosen depth
chosen_depth = 10  # This can be changed to any value you want to set as the threshold

# Step 5: Initialize a dictionary to store accuracy results for plotting
plot_data = {f'<{chosen_depth}': {df_name: [] for df_name in df_names},
             f'>={chosen_depth}': {df_name: [] for df_name in df_names}}

# Initialize a dictionary to store class-level accuracies
class_level_data = {f'<{chosen_depth}': {df_name: {class_name: [] for class_name in gene_classes} for df_name in df_names},
                    f'>={chosen_depth}': {df_name: {class_name: [] for class_name in gene_classes} for df_name in df_names}}

# Step 6: Calculate detailed accuracy for each gene in each DataFrame based on chosen depth
for df, df_name in zip(dfs, df_names):
    # Filter the DataFrame to only include the relevant genes
    relevant_genes = software_genes[df_name]
    df = df[df['Gene'].isin(relevant_genes)]

    # Calculate accuracy for each relevant gene
    for gene in relevant_genes:
        # Get the samples that are < chosen_depth and >= chosen_depth in read depth
        samples_lt_chosen = gene_reads_df[gene][gene_reads_df[gene] < chosen_depth].dropna().index
        samples_ge_chosen = gene_reads_df[gene][gene_reads_df[gene] >= chosen_depth].dropna().index

        # Filter the current DataFrame to consider only the valid samples for this gene with read depth < chosen_depth
        filtered_df_lt_chosen = df[(df['Gene'] == gene) & (df['Sample'].isin(samples_lt_chosen))]

        # Filter the current DataFrame to consider only the valid samples for this gene with read depth >= chosen_depth
        filtered_df_ge_chosen = df[(df['Gene'] == gene) & (df['Sample'].isin(samples_ge_chosen))]
        
        # Calculate accuracy for samples with read depth < chosen_depth
        if not filtered_df_lt_chosen.empty:
            total_match_lt_chosen = filtered_df_lt_chosen['Match'].sum()
            total_total_lt_chosen = filtered_df_lt_chosen['Total'].sum()
            accuracy_lt_chosen = total_match_lt_chosen / total_total_lt_chosen if total_total_lt_chosen > 0 else None
        else:
            accuracy_lt_chosen = None
            
        # Store the result in the plot_data dictionary for < chosen_depth
        plot_data[f'<{chosen_depth}'][df_name].append({'Gene': gene, 'Accuracy': accuracy_lt_chosen})

        # Calculate accuracy for samples with read depth >= chosen_depth
        if not filtered_df_ge_chosen.empty:
            total_match_ge_chosen = filtered_df_ge_chosen['Match'].sum()
            total_total_ge_chosen = filtered_df_ge_chosen['Total'].sum()
            accuracy_ge_chosen = total_match_ge_chosen / total_total_ge_chosen if total_total_ge_chosen > 0 else None
        else:
            accuracy_ge_chosen = None
            
        # Store the result in the plot_data dictionary for >= chosen_depth
        plot_data[f'>={chosen_depth}'][df_name].append({'Gene': gene, 'Accuracy': accuracy_ge_chosen})

        # Accumulate class-level data
        for class_name, class_genes in gene_classes.items():
            if gene in class_genes:
                class_level_data[f'<{chosen_depth}'][df_name][class_name].append(accuracy_lt_chosen)
                class_level_data[f'>={chosen_depth}'][df_name][class_name].append(accuracy_ge_chosen)

# Step 7: Print the gene-level results and the average accuracy for each software
for depth_category in plot_data:
    print(f"\nResults for depth category: {depth_category}")
    for df_name in df_names:
        print(f"\nSoftware: {df_name}")
        accuracies = []
        for entry in plot_data[depth_category][df_name]:
            gene = entry['Gene']
            accuracy = entry['Accuracy']
            accuracies.append(accuracy)
            if accuracy is not None:
                print(f"  Gene: {gene}, Accuracy: {accuracy:.2f}")
            else:
                print(f"  Gene: {gene}, Accuracy: No valid data")
        
        # Calculate and print the average accuracy for this software
        valid_accuracies = [acc for acc in accuracies if acc is not None]
        if valid_accuracies:
            average_accuracy = sum(valid_accuracies) / len(valid_accuracies)
            print(f"\n  Average Accuracy for {df_name}: {average_accuracy:.2f}")
        else:
            print(f"\n  Average Accuracy for {df_name}: No valid data")

# Step 8: Print the class-level accuracy results
for depth_category in class_level_data:
    print(f"\nClass-Level Results for depth category: {depth_category}")
    for df_name in df_names:
        print(f"\nSoftware: {df_name}")
        for class_name, accuracies in class_level_data[depth_category][df_name].items():
            valid_accuracies = [acc for acc in accuracies if acc is not None]
            if valid_accuracies:
                average_class_accuracy = sum(valid_accuracies) / len(valid_accuracies)
                print(f"  {class_name}: Average Accuracy = {average_class_accuracy:.2f}")
            else:
                print(f"  {class_name}: No valid data")