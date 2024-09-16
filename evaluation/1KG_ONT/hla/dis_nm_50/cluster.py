import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score

# Step 1: 读取基因距离数据
file_path = 'distance.nm.csv'
data = pd.read_csv(file_path)

# Step 2: 读取 Superpopulation 信息文件
superpop_file = 'distance.sub.nm.samples'
superpop_data = pd.read_csv(superpop_file)

# 创建样本名到 Superpopulation 的映射
sample_to_superpop = dict(zip(superpop_data['Sample'], superpop_data['Superpopulation']))

# Step 3: 处理数据
samples = pd.unique(data[['Sample1', 'Sample2']].values.ravel('K'))

# 创建样本特征矩阵
# 将每两个样本的基因距离作为特征矩阵的一行
features = []
sample_pairs = []

for i, row in data.iterrows():
    sample1 = row['Sample1']
    sample2 = row['Sample2']
    
    # 提取基因距离列（从第3列开始）
    gene_distances = row[2:].values
    features.append(gene_distances)
    sample_pairs.append((sample1, sample2))

# 将特征矩阵转换为 NumPy 数组
features = np.array(features)

# Step 4: 使用 K-Means 聚类
# 假设我们先尝试分为 3 类
kmeans = KMeans(n_clusters=3, random_state=42)
kmeans.fit(features)
cluster_labels = kmeans.labels_

# Step 5: 使用 PCA 将高维数据降维到2维
pca = PCA(n_components=2)
pca_result = pca.fit_transform(features)

# Step 6: 将 PCA 结果和聚类标签保存到 DataFrame
pca_df = pd.DataFrame(pca_result, columns=['PC1', 'PC2'])
pca_df['Cluster'] = cluster_labels
pca_df['Sample1'] = [s[0] for s in sample_pairs]
pca_df['Sample2'] = [s[1] for s in sample_pairs]
pca_df['Superpopulation1'] = pca_df['Sample1'].map(sample_to_superpop)
pca_df['Superpopulation2'] = pca_df['Sample2'].map(sample_to_superpop)

# Step 7: 可视化 PCA 降维后的聚类结果, 按照 Superpopulation 着色
plt.figure(figsize=(10, 8))

sns.scatterplot(
    x='PC1', y='PC2', 
    hue='Superpopulation1',  # 按照 Superpopulation1 着色
    style='Cluster',         # 按照聚类结果区分点的形状
    data=pca_df, 
    palette='tab10', 
    s=100, alpha=0.7
)

# 添加标题和标签
plt.title("K-Means Clustering with PCA (Colored by Superpopulation)", fontsize=16)
plt.xlabel("Principal Component 1", fontsize=12)
plt.ylabel("Principal Component 2", fontsize=12)

plt.legend(title='Superpopulation', loc='best')
plt.tight_layout()
plt.show()

# Step 8: 通过肘部法评估聚类数
sse = []
silhouette_scores = []
cluster_range = range(2, 11)  # 选择2到10个聚类数

for k in cluster_range:
    kmeans = KMeans(n_clusters=k, random_state=42)
    kmeans.fit(features)
    sse.append(kmeans.inertia_)  # SSE 评估
    silhouette_scores.append(silhouette_score(features, kmeans.labels_))  # 轮廓系数
    
# Step 9: 绘制肘部图 (Elbow Plot)
plt.figure(figsize=(10, 6))
plt.plot(cluster_range, sse, marker='o')
plt.title("Elbow Method for Optimal K", fontsize=16)
plt.xlabel("Number of Clusters", fontsize=12)
plt.ylabel("Sum of Squared Distances (SSE)", fontsize=12)
plt.show()

# Step 10: 绘制轮廓系数图 (Silhouette Score Plot)
plt.figure(figsize=(10, 6))
plt.plot(cluster_range, silhouette_scores, marker='o')
plt.title("Silhouette Score for Optimal K", fontsize=16)
plt.xlabel("Number of Clusters", fontsize=12)
plt.ylabel("Silhouette Score", fontsize=12)
plt.show()