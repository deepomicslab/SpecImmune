import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

# Step 1: 读取基因距离数据
file_path = 'distance.csv'
data = pd.read_csv(file_path)

# Step 2: 获取所有样本名
samples = pd.unique(data[['Sample1', 'Sample2']].values.ravel('K'))

# Step 3: 读取 Superpopulation 信息文件
superpop_file = 'distance.sub.samples'
superpop_data = pd.read_csv(superpop_file)

# 创建样本名到 Superpopulation 的映射
sample_to_superpop = dict(zip(superpop_data['Sample'], superpop_data['Superpopulation']))

# Step 4: 创建一个空的 N x N 距离矩阵
n = len(samples)
distance_matrix = np.zeros((n, n))

# 为样本名创建索引
sample_to_index = {sample: idx for idx, sample in enumerate(samples)}

# Step 5: 计算每对样本的欧式距离并填充矩阵
for i, row in data.iterrows():
    sample1 = row['Sample1']
    sample2 = row['Sample2']
    
    # 提取基因距离列，从 'HLA-A' 开始
    gene_distances = row[2:3].values
    
    # 计算欧式距离
    dist = np.linalg.norm(gene_distances)
    
    # 获取样本索引并在矩阵中填充
    idx1 = sample_to_index[sample1]
    idx2 = sample_to_index[sample2]
    
    distance_matrix[idx1, idx2] = dist
    distance_matrix[idx2, idx1] = dist  # 距离矩阵是对称的

# Step 6: 根据 Superpopulation 对样本排序
# 获取样本对应的 Superpopulation 标签
superpop_labels = [sample_to_superpop.get(sample, 'Unknown') for sample in samples]

# 将 Superpopulation 标签和样本名打包并按 Superpopulation 排序
sorted_samples_and_labels = sorted(zip(samples, superpop_labels), key=lambda x: x[1])
sorted_samples = [item[0] for item in sorted_samples_and_labels]
sorted_superpop_labels = [item[1] for item in sorted_samples_and_labels]

# 根据排序后的样本名重新排列距离矩阵
sorted_indices = [sample_to_index[sample] for sample in sorted_samples]
sorted_distance_matrix = distance_matrix[sorted_indices, :][:, sorted_indices]

# Step 7: 绘制热图
plt.figure(figsize=(12, 10))

# 使用 seaborn 的 heatmap 绘制距离矩阵热图
ax = sns.heatmap(
    sorted_distance_matrix, 
    xticklabels=sorted_samples, 
    yticklabels=sorted_samples, 
    cmap='coolwarm',  # 颜色映射
    cbar_kws={'label': 'Distance'},  # 色条标签
    square=True  # 正方形单元格
)

# 添加标题
plt.title("Heatmap of Sample Distance Matrix (Sorted by Superpopulation)", fontsize=16)

# 为 X 轴和 Y 轴加标签
plt.xlabel("Samples (Sorted by Superpopulation)", fontsize=12)
plt.ylabel("Samples (Sorted by Superpopulation)", fontsize=12)

# Step 8: 添加 Superpopulation 分界线
# 计算每个 Superpopulation 的边界，用于绘制分界线
prev_label = sorted_superpop_labels[0]
boundary_indices = []

for i, label in enumerate(sorted_superpop_labels):
    if label != prev_label:
        boundary_indices.append(i)
        prev_label = label

# 在热图上绘制 Superpopulation 分界线
for boundary in boundary_indices:
    ax.axhline(boundary, color='black', linewidth=2)  # 横线
    ax.axvline(boundary, color='black', linewidth=2)  # 竖线

# 显示热图
plt.tight_layout()
plt.show()

# Step 9: 保存距离矩阵为 CSV 文件
sorted_distance_df = pd.DataFrame(sorted_distance_matrix, index=sorted_samples, columns=sorted_samples)
sorted_distance_df.to_csv('sorted_sample_distance_matrix.csv')