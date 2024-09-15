import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from scipy.spatial import ConvexHull

# Step 1: 读取基因特征数据
file_path = 'distance.nm.min.ref.csv'
data = pd.read_csv(file_path)

# Step 2: 读取 Superpopulation 信息文件
superpop_file = 'distance.sub.nm.min.ref.samples'
superpop_data = pd.read_csv(superpop_file)

# 创建样本名到 Superpopulation 的映射
sample_to_superpop = dict(zip(superpop_data['Sample'], superpop_data['Superpopulation']))

# Step 3: 提取样本名和基因特征矩阵
samples = data['Sample'].values
gene_features = data.drop(columns=['Sample']).values

# Step 4: 对每个基因特征做 Z-score 标准化
scaler = StandardScaler()
gene_features_normalized = scaler.fit_transform(gene_features)

# Step 5: PCA 降维
pca = PCA(n_components=2)
pca_result = pca.fit_transform(gene_features_normalized)

# 获取解释方差比例
explained_variance = pca.explained_variance_ratio_

# Step 6: K-Means 聚类
kmeans = KMeans(n_clusters=5, random_state=42)
cluster_labels = kmeans.fit_predict(pca_result)

# Step 7: 创建 DataFrame 保存 PCA 结果和聚类标签
pca_df = pd.DataFrame(pca_result, columns=['PC1', 'PC2'])
pca_df['Sample'] = samples
pca_df['Cluster'] = cluster_labels
pca_df['Superpopulation'] = pca_df['Sample'].map(sample_to_superpop)

# Step 8: 绘制第一幅图：根据 Superpopulation 着色的 PCA 图
plt.figure(figsize=(10, 8))
sns.scatterplot(
    x='PC1', y='PC2', 
    hue='Superpopulation', 
    palette='Set1', 
    data=pca_df, 
    s=100, alpha=0.7
)
plt.title("PCA of Normalized Gene Features (Superpopulation)", fontsize=16)
plt.xlabel(f"Principal Component 1 ({explained_variance[0]*100:.2f}% Variance)", fontsize=12)
plt.ylabel(f"Principal Component 2 ({explained_variance[1]*100:.2f}% Variance)", fontsize=12)
plt.legend(title='Superpopulation', loc='best')
plt.tight_layout()
plt.show()

# Step 9: 绘制第二幅图：根据 Superpopulation 着色，并圈出每个聚类的范围
plt.figure(figsize=(10, 8))

# 绘制根据 Superpopulation 着色的散点图
sns.scatterplot(
    x='PC1', y='PC2', 
    hue='Superpopulation', 
    palette='Set1', 
    data=pca_df, 
    s=100, alpha=0.7
)

# 绘制每个簇的凸包（Convex Hull）
for i in range(kmeans.n_clusters):
    cluster_points = pca_result[cluster_labels == i]
    
    # 如果聚类中有足够的点（至少3个点才能构成凸包）
    if len(cluster_points) >= 3:
        hull = ConvexHull(cluster_points)
        hull_points = cluster_points[hull.vertices]
        plt.fill(hull_points[:, 0], hull_points[:, 1], alpha=0.2, label=f'Cluster {i}')

# 绘制聚类中心点
centroids = kmeans.cluster_centers_
plt.scatter(centroids[:, 0], centroids[:, 1], c='black', s=200, alpha=0.8, marker='X', label='Centroids')

# 添加标题和解释方差的标签
plt.title("PCA with K-Means Clustering (Superpopulation Colors)", fontsize=16)
plt.xlabel(f"Principal Component 1 ({explained_variance[0]*100:.2f}% Variance)", fontsize=12)
plt.ylabel(f"Principal Component 2 ({explained_variance[1]*100:.2f}% Variance)", fontsize=12)
plt.legend(title='Superpopulation', loc='best')
plt.tight_layout()
plt.show()