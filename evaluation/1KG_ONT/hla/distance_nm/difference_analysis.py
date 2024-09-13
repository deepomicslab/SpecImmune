import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

# Step 1: 读取基因距离数据
file_path = 'distance.nm.csv'
data = pd.read_csv(file_path)

# Step 2: 获取所有样本名
samples = pd.unique(data[['Sample1', 'Sample2']].values.ravel('K'))

# Step 3: 读取 Superpopulation 信息文件
superpop_file = 'distance.sub.nm.samples'
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
    gene_distances = row[2:].values
    
    # 计算欧式距离
    dist = np.sqrt(np.sum(gene_distances ** 2))  # 使用欧式距离公式
    
    # 获取样本索引并在矩阵中填充
    idx1 = sample_to_index[sample1]
    idx2 = sample_to_index[sample2]
    
    distance_matrix[idx1, idx2] = dist
    distance_matrix[idx2, idx1] = dist  # 距离矩阵是对称的

# Step 6: 保留所有距离：同 Superpopulation vs 不同 Superpopulation
same_superpop_distances = []
different_superpop_distances = []

# 遍历每个样本
for i, sample in enumerate(samples):
    sample_superpop = sample_to_superpop[sample]
    
    for j, other_sample in enumerate(samples):
        if i != j:  # 不考虑与自己比较的情况
            other_superpop = sample_to_superpop[other_sample]
            distance = distance_matrix[i, j]
            
            if sample_superpop == other_superpop:
                same_superpop_distances.append(distance)
            else:
                different_superpop_distances.append(distance)

# Step 7: 进行差异分析（Mann-Whitney U 检验）
same_superpop_distances = np.array(same_superpop_distances)
different_superpop_distances = np.array(different_superpop_distances)

# 使用 Mann-Whitney U 检验来比较相同 Superpopulation 和不同 Superpopulation 的距离
stat, p_value = mannwhitneyu(same_superpop_distances, different_superpop_distances, alternative='two-sided')

# 输出检验结果
print(f"Mann-Whitney U test statistic: {stat}")
print(f"P-value: {p_value}")

# Step 8: 绘制同 Superpopulation 和不同 Superpopulation 距离的分布图
plt.figure(figsize=(8, 6))

# 创建 DataFrame 来存储距离数据
distance_data = pd.DataFrame({
    'Distance': np.concatenate([same_superpop_distances, different_superpop_distances]),
    'Group': ['Same Superpopulation'] * len(same_superpop_distances) + ['Different Superpopulation'] * len(different_superpop_distances)
})

# 绘制箱线图比较两组距离
sns.boxplot(x='Group', y='Distance', data=distance_data)

# 添加标题和标签
plt.title("Comparison of Distances Within and Between Superpopulations")
plt.ylabel("Distance")
plt.xlabel("Group")

plt.tight_layout()
plt.show()