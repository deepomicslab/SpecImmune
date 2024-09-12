import pandas as pd
import matplotlib.pyplot as plt

# 定义基因分类
class1_genes = ['A', 'B', 'C', 'E', 'F', 'G']
class1_pseudogenes = ['H', 'J', 'K', 'L', 'N', 'P', 'S', 'T', 'U', 'V', 'W', 'Y']
class2_genes = ['DRA', 'DQA1', 'DQA2', 'DQB1', 'DQB2', 'DPA1', 'DPA2', 'DPB1', 'DPB2', 'DMA', 'DMB', 'DOA', 'DOB', 'DRB1', 'DRB3', 'DRB4', 'DRB5']
non_hla_genes = ['MICA', 'MICB', 'TAP1', 'TAP2']

gene_classes = {
    'Class I Genes': class1_genes,
    'Class I Pseudogenes': class1_pseudogenes,
    'Class II Genes': class2_genes,
    'Non-HLA Genes': non_hla_genes
}

# 读取 CSV 文件
csv_file = 'average_accuracies_by_class_per_depth.csv'
data = pd.read_csv(csv_file)

# 定义颜色
colors = ["#D2D2D2", "#062565", "#0098B4"]

# 计算不考虑read depth的平均准确率（按软件分开）
def calculate_overall_average_accuracy_by_software(data):
    overall_averages = {}

    tools = data['Software'].unique()

    # 对每个工具分别计算
    for tool in tools:
        overall_averages[tool] = {}
        tool_data = data[data['Software'] == tool]

        # 对每个基因类别
        for class_name, _ in gene_classes.items():
            class_data = tool_data[tool_data['Class'] == class_name]

            # 计算四个数据集的总体平均值（不考虑深度）
            average_accuracy_ge_class = class_data['Average Accuracy >= Depth'].mean()
            average_accuracy_lt_class = class_data['Average Accuracy < Depth'].mean()

            overall_averages[tool][class_name] = {
                'Average Accuracy >= Depth': average_accuracy_ge_class,
                'Average Accuracy < Depth': average_accuracy_lt_class
            }

    return overall_averages

# 计算考虑不同深度下的平均准确率（按软件分开）
def calculate_average_accuracy_by_depth_and_software(data):
    depth_averages = {}

    tools = data['Software'].unique()

    # 获取所有深度阈值
    depth_cutoffs = data['Depth Cutoff'].unique()

    # 对每个工具分别计算
    for tool in tools:
        depth_averages[tool] = {}

        tool_data = data[data['Software'] == tool]

        # 对每个深度阈值进行操作
        for depth in depth_cutoffs:
            depth_data = tool_data[tool_data['Depth Cutoff'] == depth]
            depth_averages[tool][depth] = {}

            # 对每个基因类别进行操作
            for class_name, _ in gene_classes.items():
                class_data = depth_data[depth_data['Class'] == class_name]

                # 计算每个基因类别的平均准确率
                average_accuracy_ge_class = class_data['Average Accuracy >= Depth'].mean()
                average_accuracy_lt_class = class_data['Average Accuracy < Depth'].mean()

                depth_averages[tool][depth][class_name] = {
                    'Average Accuracy >= Depth': average_accuracy_ge_class,
                    'Average Accuracy < Depth': average_accuracy_lt_class
                }

    return depth_averages

# 绘图函数 - 画出四个数据集上三个软件在不同深度阈值下的平均表现
def plot_averaged_accuracies(data):
    tools = ['HLA*LA', 'SpecHLA', 'SpecLong']

    # 遍历每类基因进行绘图
    for gene_class in gene_classes:
        # 过滤出该基因类别的数据
        class_data = data[data['Class'] == gene_class]

        # 创建绘图
        plt.figure(figsize=(4, 3))

        # 绘制每个工具的线图
        for tool, color in zip(tools, colors):
            # 过滤出每个工具的数据
            tool_data = class_data[class_data['Software'] == tool]

            # 计算四个数据集上的平均表现
            averaged_data = tool_data.groupby('Depth Cutoff')['Average Accuracy >= Depth'].mean()

            # 绘制平均表现的线图
            plt.plot(averaged_data.index, averaged_data.values, label=tool, color=color, marker='o', linewidth=3)

        # 设置标题与标签
        plt.xlabel('Depth', fontsize=9)
        plt.ylabel('Average Accuracy', fontsize=9)
        plt.ylim(0, 1.1)

        # 显示图例
        # plt.legend(title='Tool')
        # Legend outside on the right with no background
        plt.legend(title="", loc='center left', bbox_to_anchor=(1, 0.5), fontsize=6, frameon=False)

        # 显示网格
        # plt.grid(True)

        # 显示图形
        plt.tight_layout()
        # save to svg
        plt.savefig(f'average_performance_in4data_{gene_class}.svg', layout='tight', dpi=600)
        plt.show()

# 计算不考虑read depth的平均准确率（按软件分开）
overall_averages_by_software = calculate_overall_average_accuracy_by_software(data)

# 打印不考虑read depth的平均准确率
print("Overall Average Accuracies (Not Considering Depth) by Software:")
for tool, class_averages in overall_averages_by_software.items():
    print(f"\nSoftware: {tool}")
    for class_name, averages in class_averages.items():
        print(f"  {class_name}:")
        print(f"    Average Accuracy >= Depth: {averages['Average Accuracy >= Depth']:.4f}")
        print(f"    Average Accuracy < Depth: {averages['Average Accuracy < Depth']:.4f}")

# 计算考虑不同深度下的平均准确率（按软件分开）
depth_averages_by_software = calculate_average_accuracy_by_depth_and_software(data)

# 打印考虑不同深度下的平均准确率
print("\nAverage Accuracies by Depth by Software:")
for tool, depth_averages in depth_averages_by_software.items():
    print(f"\nSoftware: {tool}")
    for depth, class_averages in depth_averages.items():
        print(f"  Depth Cutoff: {depth}")
        for class_name, averages in class_averages.items():
            print(f"    {class_name}:")
            print(f"      Average Accuracy >= Depth: {averages['Average Accuracy >= Depth']:.4f}")
            print(f"      Average Accuracy < Depth: {averages['Average Accuracy < Depth']:.4f}")

# 绘制图形，显示四个数据集上三个软件在不同深度阈值下的平均表现
plot_averaged_accuracies(data)