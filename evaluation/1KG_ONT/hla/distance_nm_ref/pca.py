import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler,  MinMaxScaler, MaxAbsScaler, RobustScaler, QuantileTransformer, PowerTransformer
from sklearn.decomposition import PCA  # Import PCA
from sklearn.feature_selection import VarianceThreshold


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
gene_features = data.drop(columns=['Sample'])

# Remove the 'HLA-DQA2' feature
if 'HLA-DRB1' in gene_features.columns:
    gene_features = gene_features.drop(columns=['HLA-DRB1'])
gene_features = gene_features.drop(columns=['HLA-A'])
gene_features = gene_features.drop(columns=['HLA-B'])
gene_features = gene_features.drop(columns=['HLA-C'])
gene_features = gene_features.drop(columns=['HLA-DQA1'])
gene_features = gene_features.drop(columns=['HLA-DQA2'])

gene_features = gene_features.drop(columns=['HLA-DPB2'])
gene_features = gene_features.drop(columns=['HLA-DQB1'])
# gene_features = gene_features.drop(columns=['HLA-DPA1'])
# gene_features = gene_features.drop(columns=['HLA-DPB1'])



print(gene_features)

# Step 5: Standardize the gene features using Z-score normalization
# scaler = StandardScaler()
# gene_features_normalized = scaler.fit_transform(gene_features)
# gene_features_log_transformed = np.log1p(gene_features_high_variance)
# Initialize MinMaxScaler
gene_features_log_transformed = np.log1p(gene_features)

# Optionally, follow this with Min-Max scaling or StandardScaler
scaler = RobustScaler()  # or MinMaxScaler()
gene_features_normalized = scaler.fit_transform(gene_features)



import umap

# Apply UMAP to the normalized data
umap_model = umap.UMAP(n_components=2, random_state=42)
umap_result = umap_model.fit_transform(gene_features_normalized)

# You can plot the results using the same method as before

# You can plot the results using the same method as before

# print(gene_features_normalized)
# Initialize PowerTransformer (Yeo-Johnson works with positive/negative data)
# scaler = PowerTransformer(method='yeo-johnson')

# # Apply Power transformation
# gene_features_normalized = scaler.fit_transform(gene_features_log_transformed)
# Step 6: PCA dimensionality reduction to 2 components
# pca = PCA(n_components=2)
# pca_result = pca.fit_transform(gene_features)

# Step 7: Create DataFrame with PCA results
pca_df = pd.DataFrame({
    'PCA1': umap_result[:, 0],
    'PCA2': umap_result[:, 1],
    'Sample': samples
})

# Map Superpopulation labels to the DataFrame
pca_df['Superpopulation'] = pca_df['Sample'].map(sample_to_superpop)
pca_df = pca_df.dropna(subset=['Superpopulation'])

# Step 8: Plot PCA results colored by Superpopulation
plt.figure(figsize=(10, 8))
sns.scatterplot(
    x='PCA1', y='PCA2',
    hue='Superpopulation',
    data=pca_df,
    palette='Set1',
    s=100,
    alpha=0.7
)
plt.title("PCA of Normalized Gene Features by Superpopulation", fontsize=16)
plt.xlabel("PCA Component 1", fontsize=12)
plt.ylabel("PCA Component 2", fontsize=12)
plt.legend(title='Superpopulation', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.show()