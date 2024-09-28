import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# 处理的基因数目
genes_count = {
    'HLA*LA': 17,
    'SpecHLA': 8,
    'SpecLong': 39
}

# 读取三个 CSV 文件，分别代表三个软件的性能数据
software1 = pd.read_csv('hlala.time.csv')
software2 = pd.read_csv('spechla.hla.time.csv')
software3 = pd.read_csv('speclong.hla.time.csv')

# 添加一列 'software' 用于区分不同软件
software1['software'] = 'HLA*LA'
software2['software'] = 'SpecHLA'
software3['software'] = 'SpecLong'

# 合并三个 DataFrame
data = pd.concat([software1, software2, software3])

# 归一化处理时间和内存
data['normalized_time'] = data.apply(lambda row: row['time'] / genes_count[row['software']], axis=1)
data['normalized_mem'] = data.apply(lambda row: row['mem'] / genes_count[row['software']], axis=1)
data['normalized_wall_clock_time'] = data.apply(lambda row: row['wall_clock_time'] / genes_count[row['software']], axis=1)

# 计算每个软件的平均时间、墙上时钟时间和内存使用
averages = data.groupby('software').agg({
    'time': 'mean',
    'wall_clock_time': 'mean',
    'mem': 'mean',
    'normalized_time': 'mean',
    'normalized_wall_clock_time': 'mean',
    'normalized_mem': 'mean'
}).reset_index()

# 打印平均值
print("Average values per software:")
print(averages)

# 设置 Seaborn 风格
sns.set(style="whitegrid")

# 创建一个包含 3 行 1 列的图表布局
fig, axes = plt.subplots(3, 1, figsize=(20, 12))

# 1. 可视化归一化后的墙上时钟时间 (normalized_time) 的比较
barplot1 = sns.barplot(x="sample", y="normalized_wall_clock_time", hue="software", data=data, ax=axes[0])
axes[0].set_title('Normalized Wall Clock Time Comparison (per gene)')
axes[0].set_xlabel('Sample')
axes[0].set_ylabel('Normalized Time (hours per gene)')
# rotate x-axis labels
axes[0].tick_params(axis='x', rotation=45)

# 2. 可视化归一化后的内存使用 (normalized_mem) 的比较
barplot2 = sns.barplot(x="sample", y="mem", hue="software", data=data, ax=axes[1])
axes[1].set_title('Memory Usage Comparison')
axes[1].set_xlabel('Sample')
axes[1].set_ylabel('Memory Usage (GB)')
# rotate x-axis labels
axes[1].tick_params(axis='x', rotation=45)

# 3. 可视化归一化后的时间 (normalized_time) 的比较
barplot3 = sns.barplot(x="sample", y="normalized_time", hue="software", data=data, ax=axes[2])
axes[2].set_title('Normalized Time Comparison (per gene)')
axes[2].set_xlabel('Sample')
axes[2].set_ylabel('Normalized Time (hours per gene)')
# rotate x-axis labels
axes[2].tick_params(axis='x', rotation=45)

# Remove individual legends in subplots
axes[0].legend_.remove()
axes[1].legend_.remove()
axes[2].legend_.remove()

# Extract handles and labels from one of the barplots to keep the colors consistent
handles, labels = axes[0].get_legend_handles_labels()

# Add a single legend below the entire figure
fig.legend(handles=handles, labels=labels, loc='lower center', ncol=3, bbox_to_anchor=(0.5, -0.05))

# 调整布局并显示图表
plt.tight_layout()
plt.savefig('compare_hla_resource.pdf', format='pdf', dpi=600, bbox_inches='tight')
plt.show()