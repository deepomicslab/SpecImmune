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



# Use the specified colors
colors = ["#D2D2D2", "#062565", "#0098B4"]
# Step 7: Calculate and plot average accuracy across all genes and samples for each depth cutoff

# Initialize a dictionary to store the average accuracy for each DataFrame at each cutoff
average_accuracies = {df_name: [] for df_name in df_names}

for cutoff in read_cutoffs:
    for df_name in df_names:
        # Extract accuracies for all genes at this cutoff
        accuracies = [item['Accuracy'] for item in plot_data[cutoff][df_name] if item['Accuracy'] is not None]
        
        # Calculate the average accuracy (if there are any accuracies to average)
        avg_accuracy = sum(accuracies) / len(accuracies) if accuracies else None
        
        # Store the average accuracy
        average_accuracies[df_name].append(avg_accuracy)

# Create the plot for average accuracy vs. depth cutoff
plt.figure(figsize=(4, 3))  # Adjust the figure size

for df_name, color in zip(df_names, colors):
    plt.plot(read_cutoffs, average_accuracies[df_name], label=df_name, color=color, marker='o', linewidth=3)

# Add plot labels and titles
plt.xlabel("Depth", fontsize=9)
plt.ylabel("Average Accuracy", fontsize=9)
plt.ylim(0, 1.1)  # Assuming accuracy is between 0 and 1
# plt.xticks(read_cutoffs, fontsize=9)
# plt.yticks(fontsize=9)

# Legend outside on the right with no background
plt.legend(title="", loc='center left', bbox_to_anchor=(1, 0.5), fontsize=6, frameon=False)

plt.tight_layout(pad=2)  # Adjust padding for better spacing
plt.savefig("overview_average_depth_accuracy_hla.svg", format="svg", bbox_inches='tight', dpi=600)  # Save the overview plot as an SVG
plt.show()  # Display the plot

