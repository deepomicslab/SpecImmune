import pandas as pd
import matplotlib.pyplot as plt
import argparse

def plot_gtf_first_n_genes(gtf_file, n=10):
    # Read the GTF file
    columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
    gtf = pd.read_csv(gtf_file, sep='\t', comment='#', names=columns)
    
    # Extract gene names and positions
    gtf['gene_id'] = gtf['attribute'].str.extract('gene_id "([^"]+)"')
    gtf['gene_name'] = gtf['attribute'].str.extract('gene_name "([^"]+)"')
    genes = gtf[['gene_id', 'gene_name', 'seqname', 'start', 'end']].drop_duplicates()
    
    # Select the first N genes
    first_n_genes = genes.head(n)
    
    if first_n_genes.empty:
        print(f"No genes found in the GTF file.")
        return
    
    # Determine the start and end positions dynamically
    start = first_n_genes['start'].min()
    end = first_n_genes['end'].max()
    chromosome = first_n_genes['seqname'].iloc[0]
    
    # Plot the region
    plt.figure(figsize=(10, 3))
    y_positions = range(1, len(first_n_genes) + 1)  # Assign unique y-positions to each gene
    
    for y, (_, gene) in zip(y_positions, first_n_genes.iterrows()):
        plt.plot([gene['start'], gene['end']], [y, y], label=f"{gene['gene_name']} ({gene['gene_id']})")
    
    plt.xlabel('Genomic Position')
    plt.ylabel('Genes')
    plt.title(f'First {n} Genes: Chromosome {chromosome} ({start}-{end})')
    plt.yticks(y_positions, first_n_genes['gene_name'])
    plt.legend(loc='upper right', bbox_to_anchor=(1.3, 1))
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot the first N genes from a GTF file.')
    parser.add_argument('gtf_file', type=str, help='Path to the GTF file')
    parser.add_argument('--num_genes', type=int, default=10, help='Number of genes to plot')
    
    args = parser.parse_args()
    
    plot_gtf_first_n_genes(args.gtf_file, args.num_genes)