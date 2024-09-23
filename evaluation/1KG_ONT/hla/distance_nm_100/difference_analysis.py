import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns

# Step 1: 读取基因距离数据
file_path = 'distance.nm.min100.csv'  # 替换为你的文件路径
data = pd.read_csv(file_path)

# Step 2: 读取 Superpopulation 信息文件
superpop_file = 'distance.sub.nm.min100.samples'  # 替换为你的文件路径
superpop_data = pd.read_csv(superpop_file)

# 创建样本名到 Superpopulation 的映射
sample_to_superpop = dict(zip(superpop_data['Sample'], superpop_data['Superpopulation']))

# Step 3: 创建两组：相同 Superpopulation 和 不同 Superpopulation
same_superpop_distances = []
diff_superpop_distances = []

for index, row in data.iterrows():
    sample1 = row['Sample1']
    sample2 = row['Sample2']
    
    # 提取两样本的 Superpopulation
    superpop1 = sample_to_superpop.get(sample1)
    superpop2 = sample_to_superpop.get(sample2)
    
    # 提取基因距离（从第3列开始）
    gene_distance = np.median(row[11:12]) # 假设基因距离列是从第3列开始的，取均值
    
    # 将距离归类到相同或不同 Superpopulation 组
    if superpop1 == superpop2:
        same_superpop_distances.append(gene_distance)
    else:
        diff_superpop_distances.append(gene_distance)

# 转换为 NumPy 数组
same_superpop_distances = np.array(same_superpop_distances)
diff_superpop_distances = np.array(diff_superpop_distances)

# Step 4: 可视化两个组的基因距离分布
plt.figure(figsize=(10, 6))

sns.histplot(same_superpop_distances, color='blue', label='Same Superpopulation', kde=True, stat="density", bins=30)
sns.histplot(diff_superpop_distances, color='red', label='Different Superpopulation', kde=True, stat="density", bins=30)

plt.title("Distribution of Genetic Distances (Same vs Different Superpopulation)", fontsize=16)
plt.xlabel("Genetic Distance", fontsize=12)
plt.ylabel("Density", fontsize=12)
plt.legend(loc='upper right')
plt.tight_layout()
# save
# plt.savefig('distance_distribution.svg', format='svg', bbox_inches='tight', dpi=600)
plt.savefig('distance_distribution_DQA2.pdf', format='pdf', bbox_inches='tight', dpi=600)

plt.show()

# Step 5: 使用 Mann-Whitney U 检验比较两组之间的差异
u_statistic, p_value_mannwhitney = stats.mannwhitneyu(same_superpop_distances, diff_superpop_distances, alternative='two-sided')

# Step 5b: 使用 t-test 比较两组之间的差异（如果数据是正态分布）
t_statistic, p_value_ttest = stats.ttest_ind(same_superpop_distances, diff_superpop_distances)

# 输出 p-value 结果
print(f"Mann-Whitney U test p-value: {p_value_mannwhitney:.6f}")
print(f"T-test p-value: {p_value_ttest:.6f}")

# 可选：输出统计量
print(f"Mann-Whitney U test statistic: {u_statistic}")
print(f"T-test statistic: {t_statistic}")