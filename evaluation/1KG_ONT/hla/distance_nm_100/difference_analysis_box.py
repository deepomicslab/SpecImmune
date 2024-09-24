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


required_genes = ['HLA-A', 'HLA-B', 'HLA-C', 'HLA-DPA1', 'HLA-DPB1', 'HLA-DQB1']
missing_genes = [gene for gene in required_genes if gene not in data.columns]
# Step 3: 创建两组：相同 Superpopulation 和 不同 Superpopulation
same_superpop_distances = []
diff_superpop_distances = []

for index, row in data.iterrows():
    sample1 = row['Sample1']
    sample2 = row['Sample2']
    
    # 提取两样本的 Superpopulation
    superpop1 = sample_to_superpop.get(sample1)
    superpop2 = sample_to_superpop.get(sample2)
    
    # 如果 Superpopulation 信息缺失，跳过该行
    if superpop1 is None or superpop2 is None:
        continue
    
    # 提取基因距离（从第3列开始）
    gene_distances = row[required_genes].values  # 或者使用 .loc[row.index, required_genes]
    gene_distance = np.median(gene_distances)  # 取中位数
    # gene_distance = np.median(row[2:])  # 假设基因距离列是从第3列开始的，取中位数
    
    # 将距离归类到相同或不同 Superpopulation 组
    if superpop1 == superpop2:
        same_superpop_distances.append(gene_distance)
    else:
        diff_superpop_distances.append(gene_distance)

# 转换为 NumPy 数组
same_superpop_distances = np.array(same_superpop_distances)
diff_superpop_distances = np.array(diff_superpop_distances)

# Step 4: 可视化两个组的基因距离分布（箱线图）

# 创建一个 DataFrame 用于绘制箱线图
df = pd.DataFrame({
    'Distance': np.concatenate([same_superpop_distances, diff_superpop_distances]),
    'Group': ['Intra-superpopulation'] * len(same_superpop_distances) + ['Inter-superpopulation'] * len(diff_superpop_distances)
})

# out to a file
df.to_csv('distance_comparison_inter_intra.csv', index=False)

# 创建箱线图，去除离群点
plt.figure(figsize=(5, 5))

# 自定义调色板
custom_palette = ['#7F7FFF', '#FF7F7F']  # 蓝色和红色

# 绘制箱线图，调整风格
ax = sns.boxplot(
    x='Group', 
    y='Distance', 
    data=df, 
    width=0.3,
    showfliers=False,  # 不显示离群点
    palette=custom_palette,  # 自定义调色板
    linewidth=1.5,     # 设置箱子的边框厚度
    boxprops=dict(edgecolor='black'),  # 自定义箱子的边框颜色
    whiskerprops=dict(color='black', linewidth=1.5),  # 须线风格
    capprops=dict(color='none'),  # 去掉顶端线
    medianprops=dict(color='black', linewidth=2.5)  # 中位数线的风格
)

# Step 5: 使用 t 检验比较两组之间的差异
t_statistic, p_value_ttest = stats.ttest_ind(
    same_superpop_distances, 
    diff_superpop_distances, 
    equal_var=False  # 使用 Welch's t-test，不假定方差相等
)

print(f"T-statistic: {t_statistic}")
print(f"Degrees of freedom: {len(same_superpop_distances) + len(diff_superpop_distances) - 2}")
# 获取 y 轴的上下限，以便确定 p 值显示的位置
ymin, ymax = ax.get_ylim()

# 在图表上方显示 p 值（科学计数法）
p_value_text = f"T-test p-value: {p_value_ttest:.2e}"
# plt.text(0.5, ymax * 0.95, p_value_text, ha='center', va='bottom', fontsize=12, color='black')

# 设置标题和标签
plt.title("Genetic Distance Comparison\n(Inter vs Intra Superpopulation)", fontsize=12, weight='bold')
plt.xlabel("Group", fontsize=12, weight='bold')
plt.ylabel("Genetic Distance", fontsize=12, weight='bold')

# 优化布局
plt.tight_layout()

# 保存图表为 PDF 文件
plt.savefig('genetic_distance_comparison.pdf', format='pdf', dpi=600, bbox_inches='tight')

# 显示图表
plt.show()

# 输出 p-value 结果（科学计数法）
print(f"T-test p-value: {p_value_ttest:.2e}")