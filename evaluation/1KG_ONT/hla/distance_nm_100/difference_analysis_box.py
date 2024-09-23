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
    
    # 如果 Superpopulation 信息缺失，跳过该行
    if superpop1 is None or superpop2 is None:
        continue
    
    # 提取基因距离（从第3列开始）
    gene_distance = np.median(row[11:12])  # 假设基因距离列是从第3列开始的，取中位数
    
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
    'Group': ['Same Superpopulation'] * len(same_superpop_distances) + ['Different Superpopulation'] * len(diff_superpop_distances)
})

# 创建箱线图，去除离群点
plt.figure(figsize=(8, 6))

# 这里我们使用 `sns.histplot` 中的蓝色和红色
custom_palette = ['#7F7FFF', '#FF7F7F']  # 蓝色 (blue) 和 红色 (red)

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

# Step 5: 使用 Mann-Whitney U 检验比较两组之间的差异
u_statistic, p_value_mannwhitney = stats.mannwhitneyu(same_superpop_distances, diff_superpop_distances, alternative='two-sided')

# 获取 y 轴的上下限，以便确定 p 值显示的位置
ymin, ymax = ax.get_ylim()

# 在水平线的上方显示 p-value（去掉横线）
p_value_text = f"Mann-Whitney U test p-value: {p_value_mannwhitney:.6f}"
plt.text(0.5, ymax * 0.95, p_value_text, ha='center', va='bottom', fontsize=12, color='black')

# 设置标题和标签
plt.title("Genetic Distance Comparison\n(Same vs Different Superpopulation)", fontsize=16, weight='bold')
plt.xlabel("Group", fontsize=12, weight='bold')
plt.ylabel("Genetic Distance", fontsize=12, weight='bold')

# 优化布局
plt.tight_layout()

# save the plot to svg
# plt.savefig('genetic_distance_comparison.svg', format='svg', dpi=600, box_inches='tight')
plt.savefig('genetic_distance_comparison_DQA2.pdf', format='pdf', dpi=600, box_inches='tight')

# 显示图表
plt.show()

# 输出 p-value 结果
print(f"Mann-Whitney U test p-value: {p_value_mannwhitney:.6f}")