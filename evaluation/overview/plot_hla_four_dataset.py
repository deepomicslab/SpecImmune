import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# 定义基因分类
class1_genes = ['A', 'B', 'C', 'E', 'F', 'G']
class1_pseudogenes = ['H', 'J', 'K', 'L', 'N', 'P', 'S', 'T', 'U', 'V', 'W', 'Y']
class1_pseudogenes = ['H', 'J', 'K', 'L', 'N', 'P', 'S', 'T', 'U', 'V', 'W']

class2_genes = ['DRA', 'DQA1', 'DQA2', 'DQB1', 'DQB2', 'DPA1', 'DPA2', 'DPB1', 'DPB2', 'DMA', 'DMB', 'DOA', 'DOB', 'DRB1', 'DRB3', 'DRB4', 'DRB5']
class2_genes = ['DRA', 'DQA1', 'DQA2', 'DQB1', 'DPA1', 'DPA2', 'DPB1', 'DPB2', 'DMA', 'DMB', 'DOA', 'DOB', 'DRB1', 'DRB3']

non_hla_genes = ['MICA', 'MICB', 'TAP1', 'TAP2']

# 合并所有基因分类
gene_order = class1_genes + class1_pseudogenes + class2_genes + non_hla_genes
colors = ["#D2D2D2", "#062565", "#0098B4"]

# 计算准确率
def calculate_accuracy(df):
    df['Accuracy'] = df['Match'] / df['Total']
    return df

# 从文件中读取数据
def load_data_from_files(dataset_dirs):
    all_data = {}
    for dataset, data_dir in dataset_dirs.items():
        all_data[dataset] = {}
        for software in ['hlala', 'spechla', 'speclong']:
            file_path = os.path.join(data_dir, f'{software}.match.csv')
            if os.path.exists(file_path):
                df = pd.read_csv(file_path)
                all_data[dataset][software] = df
            else:
                print(f"File {file_path} not found.")
    return all_data

def plot_accuracy(all_data):
    for dataset, software_data in all_data.items():
        combined_data = []

        for software, df in software_data.items():
            # 计算准确率
            df = calculate_accuracy(df)
            
            # 计算每个基因的平均准确率（跨多个样本）
            gene_accuracy = df.groupby('Gene')['Accuracy'].mean().reset_index()
            gene_accuracy['Software'] = software  # 添加软件名称列
            
            # 按照基因分类排序
            gene_accuracy['Gene'] = pd.Categorical(gene_accuracy['Gene'], categories=gene_order, ordered=True)
            gene_accuracy = gene_accuracy.sort_values('Gene')
            
            # 将结果添加到 combined_data 列表
            combined_data.append(gene_accuracy)
        
        # 合并三个软件的数据
        combined_df = pd.concat(combined_data)
        
        # 绘制柱状图
        plt.figure(figsize=(12, 3))  # 调整图像大小，减少高度
        sns.barplot(x='Gene', y='Accuracy', hue='Software', data=combined_df, palette=colors, ci=None, dodge=True)
        
        # plt.title(f'{dataset}: Gene Accuracy Comparison Across Software')
        plt.xticks(rotation=45)
        plt.xlabel('Locus',fontsize=14)
        plt.ylabel('Accuracy',fontsize=14)
        # set label font size
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.ylim(0, 1.1)
        plt.legend( 
            loc='center left', 
            bbox_to_anchor=(1, 0.5), 
            fontsize=14, 
            frameon=False
        )

        plt.tight_layout()
        
        # 保存图像为SVG格式
        svg_file_path = f'{dataset}_gene_accuracy_comparison.svg'
        plt.savefig(svg_file_path, format='svg', bbox_inches='tight', dpi=600)

        # 显示图像
        plt.show()

# 指定每个数据集的目录
dataset_dirs = {
    'hgsvc2_clr': '../HGSVC2/clr',  # 修改为实际路径
    'hgsvc2_hifi': '../HGSVC2/hifi',  # 修改为实际路径
    'hprc_hifi': '../HPRC/hla/hifi',  # 修改为实际路径
    'hprc_ont': '../HPRC/hla/ont'  # 修改为实际路径
}

# 读取文件并加载数据
all_data = load_data_from_files(dataset_dirs)

# 绘制准确率柱状图并保存为SVG
plot_accuracy(all_data)