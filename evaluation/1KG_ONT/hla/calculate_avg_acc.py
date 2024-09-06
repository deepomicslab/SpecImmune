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

# Step 2: Read the gene read counts file into a DataFrame
gene_reads_df = pd.read_csv('read_depth.csv', index_col=0)
gene_reads_df.columns = gene_reads_df.columns.str.replace('HLA-', '', regex=False)

# Filter gene_reads_df to only include genes present in any of the DataFrames
genes_in_dfs = pd.concat([df['Gene'] for df in dfs]).unique()
gene_reads_df = gene_reads_df[genes_in_dfs]

# Sort genes alphabetically
genes_in_dfs.sort()

# Step 3: Set the chosen depth
chosen_depth = 15  # This can be changed to any value you want to set as the threshold

# Step 4: Initialize a dictionary to store accuracy results for plotting
plot_data = {f'<{chosen_depth}': {df_name: [] for df_name in df_names},
             f'>={chosen_depth}': {df_name: [] for df_name in df_names}}

# Step 5: Calculate detailed accuracy for each gene in each DataFrame based on chosen depth
for df, df_name in zip(dfs, df_names):
    # Calculate accuracy for each gene
    for gene in genes_in_dfs:
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

# Step 6: Print the gene-level results and the average accuracy for each software
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