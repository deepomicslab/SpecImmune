import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

# Step 1: 读取基因距离数据
file_path = 'distance.nm.min.csv'
data = pd.read_csv(file_path)

# Step 2: 读取 Superpopulation 信息文件
superpop_file = 'distance.sub.nm.min.samples'
superpop_data = pd.read_csv(superpop_file)

# 创建样本名到 Superpopulation 的映射
sample_to_superpop = dict(zip(superpop_data['Sample'], superpop_data['Superpopulation']))

# Step 3: 处理数据
# 假设每个样本在 Sample1 和 Sample2 中都出现，我们提取所有样本名
samples = pd.unique(data[['Sample1', 'Sample2']].values.ravel('K'))

# 创建样本特征矩阵
# 将每两个样本的基因距离的中位数作为特征
features = []
sample_pairs = []

for i, row in data.iterrows():
    sample1 = row['Sample1']
    sample2 = row['Sample2']
    
    # 提取基因距离列（从第3列开始）
    gene_distances = row[2:].values
    
    # 计算基因距离的中位数
    median_distance = np.median(gene_distances)
    
    # 将中位数作为特征
    features.append([median_distance])  # 注意：这里是一个单一元素的列表
    sample_pairs.append((sample1, sample2))

# 将特征矩阵转换为 NumPy 数组
features = np.array(features)

# Step 4: PCA 降维
# 由于我们只有一个特征（中位数距离），PCA 可能不是必要的，但为了展示，我们将它保留
pca = PCA(n_components=1)  # 降到2维
pca_result = pca.fit_transform(features)

# Step 5: 创建 DataFrame 保存 PCA 结果
pca_df = pd.DataFrame(pca_result, columns=['PC1', 'PC2'])
pca_df['Sample1'] = [s[0] for s in sample_pairs]
pca_df['Sample2'] = [s[1] for s in sample_pairs]

# 将 Superpopulation 标签添加到每个样本
pca_df['Superpopulation1'] = pca_df['Sample1'].map(sample_to_superpop)
pca_df['Superpopulation2'] = pca_df['Sample2'].map(sample_to_superpop)

# Step 6: 可视化 PCA 结果
plt.figure(figsize=(10, 8))

# 我们只绘制 Sample1 和其 Superpopulation 的 PCA 结果
sns.scatterplot(
    x='PC1', y='PC2', 
    hue='Superpopulation1', 
    data=pca_df, 
    palette='Set1', 
    s=100, alpha=0.7
)

# 添加标题和标签
plt.title("PCA of Gene Distances Based on Sample Pairs (Using Median Distance)", fontsize=16)
plt.xlabel("Principal Component 1", fontsize=12)
plt.ylabel("Principal Component 2", fontsize=12)

plt.legend(title='Superpopulation', loc='best')
plt.tight_layout()
plt.show()