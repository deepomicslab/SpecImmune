from sklearn.manifold import TSNE

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from scipy.spatial.distance import pdist
import seaborn as sns

# Step 1: Read gene feature data
file_path = 'distance.nm.min.ref.csv'
data = pd.read_csv(file_path)

# Step 2: Read Superpopulation information file
superpop_file = 'distance.sub.nm.min.ref.samples'
superpop_data = pd.read_csv(superpop_file)

# Step 3: Create sample to superpopulation mapping
sample_to_superpop = dict(zip(superpop_data['Sample'], superpop_data['Superpopulation']))

# Step 4: Extract sample names and gene feature matrix
samples = data['Sample'].values
gene_features = data.drop(columns=['Sample']).values

# Step 5: Standardize the gene features using Z-score normalization
scaler = StandardScaler()
gene_features_normalized = scaler.fit_transform(gene_features)

# Step 6: Try different perplexity values
perplexities = [5, 30, 50, 100]

for perplexity in perplexities:
    tsne = TSNE(n_components=2, perplexity=perplexity, random_state=42)
    tsne_result = tsne.fit_transform(gene_features_normalized)
    
    tsne_df = pd.DataFrame({
        'Dim1': tsne_result[:, 0],
        'Dim2': tsne_result[:, 1],
        'Sample': samples
    })
    tsne_df['Superpopulation'] = tsne_df['Sample'].map(sample_to_superpop)
    tsne_df = tsne_df.dropna(subset=['Superpopulation'])
    
    plt.figure(figsize=(10, 8))
    sns.scatterplot(
        x='Dim1', y='Dim2',
        hue='Superpopulation',
        data=tsne_df,
        palette='Set1',
        s=100,
        alpha=0.7
    )
    plt.title(f"t-SNE with Perplexity={perplexity}", fontsize=16)
    plt.xlabel("t-SNE Dimension 1", fontsize=12)
    plt.ylabel("t-SNE Dimension 2", fontsize=12)
    plt.legend(title='Superpopulation', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.show()