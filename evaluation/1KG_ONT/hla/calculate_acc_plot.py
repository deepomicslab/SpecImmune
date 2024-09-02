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
    print(plot_data[cutoff])


# exit()  # Remove this line to continue with the rest of the code
# Step 5: Plot the first set of results by cutoff and save to PDF
plt.figure(figsize=(10, 12))  # Adjusted figure size with smaller height for subplots

# Use the specified colors
colors = ["#FF5959", "#FFAD5A",  "#4F9DA6"]
colors = ["#D2D2D2", "#062565", "#0098B4"]

for i, cutoff in enumerate(read_cutoffs, 1):
    plt.subplot(5, 2, i)  # Create subplots in a 5x2 grid
    
    for j, (df_name, color) in enumerate(zip(df_names, colors)):
        gene_accuracies = plot_data[cutoff][df_name]
        genes = [item['Gene'] for item in gene_accuracies]
        accuracies = [item['Accuracy'] for item in gene_accuracies]
        
        sns.lineplot(x=genes, y=accuracies, label=df_name, color=color, marker='o', linewidth=3)
    
    plt.title(f"Depth Cutoff = {cutoff}", fontsize=14)
    plt.xlabel("Locus", fontsize=14)
    plt.ylabel("Accuracy", fontsize=14)
    plt.xticks(rotation=45, fontsize=14)
    plt.yticks(fontsize=14)
    plt.ylim(0, 1.1)  # Assuming accuracy is between 0 and 1

    
    # Legend outside on the right with no background
    plt.legend(title="", loc='center left', bbox_to_anchor=(1, 0.5), fontsize=14, frameon=False)

plt.tight_layout(pad=1)  # Adjust padding for better spacing
plt.savefig("cutoff_accuracy_plots.svg", format="svg", bbox_inches='tight', dpi=600)  # Save the first plot as a PDF
# plt.show()

# Step 6: Plot the second set of results by gene, each gene in a row, and save to PDF
plt.figure(figsize=(6, len(genes_in_dfs) * 2))  # Adjusted figure size with smaller width

for i, gene in enumerate(genes_in_dfs, 1):
    plt.subplot(len(genes_in_dfs), 1, i)  # Each gene gets its own row
    
    for j, (df_name, color) in enumerate(zip(df_names, colors)):
        accuracies = [plot_data[cutoff][df_name][i-1]['Accuracy'] if i-1 < len(plot_data[cutoff][df_name]) else None for cutoff in read_cutoffs]
        
        sns.lineplot(x=read_cutoffs, y=accuracies, label=df_name, color=color, marker='o')
    
    plt.title(f"{gene}", fontsize=12)
    plt.xlabel("Depth Cutoff", fontsize=14)
    plt.ylabel("Accuracy", fontsize=10)
    plt.xticks(rotation=45, fontsize=8)
    plt.yticks(fontsize=14)
    plt.ylim(0, 1.1)  # Assuming accuracy is between 0 and 1
    
    # Legend outside on the right with no background
    plt.legend(title="", loc='center left', bbox_to_anchor=(1, 0.5), fontsize=8, frameon=False)

plt.tight_layout(pad=1)  # Adjust padding for better spacing
plt.savefig("gene_accuracy_by_cutoff_plots.svg", format="svg", bbox_inches='tight', dpi=600)  # Save the second plot as a PDF
# plt.show()