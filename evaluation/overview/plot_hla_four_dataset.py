import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import numpy as np

# 定义基因分类
class1_genes = ['A', 'B', 'C', 'E', 'F', 'G']
class1_pseudogenes = ['H', 'J', 'K', 'L', 'N', 'P', 'S', 'T', 'U', 'V', 'W', 'Y']

class2_genes = ['DRA', 'DQA1', 'DQA2', 'DQB1', 'DQB2', 'DPA1', 'DPA2', 'DPB1', 'DPB2', 'DMA', 'DMB', 'DOA', 'DOB', 'DRB1', 'DRB3', 'DRB4', 'DRB5']
non_hla_genes = ['MICA', 'MICB', 'TAP1', 'TAP2']

# 合并所有基因分类
gene_order = class1_genes + class1_pseudogenes + class2_genes + non_hla_genes

KGP_genes = ['A', 'B', 'C', 'DQB1', 'DRB1']
colors = ["#D2D2D2", "#062565", "#0098B4"]

# 软件名称映射字典
software_name_map = {
    'hlala': 'HLA*LA',
    'spechla': 'SpecHLA',
    'speclong': 'SpecLong'
}

# 计算准确率
def calculate_accuracy(df):
    df['Accuracy'] = df['Match'] / df['Total']
    return df

# 从文件中读取数据，包括深度信息
def load_data_from_files(dataset_dirs):
    all_data = {}
    for dataset, data_dir in dataset_dirs.items():
        all_data[dataset] = {'software_data': {}, 'depth_data': None}
        for software in ['hlala', 'spechla', 'speclong']:
            file_path = os.path.join(data_dir, f'{software}.match.csv')
            if os.path.exists(file_path):
                df = pd.read_csv(file_path)
                all_data[dataset]['software_data'][software] = df
            else:
                print(f"File {file_path} not found.")
        
        # Load depth data for the dataset
        depth_file_path = os.path.join(data_dir, 'read_depth.csv')
        if os.path.exists(depth_file_path):
            depth_df = pd.read_csv(depth_file_path)
            all_data[dataset]['depth_data'] = depth_df
        else:
            print(f"Depth file {depth_file_path} not found.")
    
    return all_data

import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# 绘图函数，包含准确率和带有雨云图的深度图
def plot_accuracy(all_data):
    for dataset, data in all_data.items():
        software_data = data['software_data']
        depth_data = data['depth_data']

        # 如果数据集是 '1kg_ont'，只聚焦于 KGP_genes，并缩小图的宽度
        if dataset == '1kg_ont':
            gene_filter = KGP_genes
            fig_width = 4  # 缩小宽度
        else:
            gene_filter = gene_order
            fig_width = 12  # 默认宽度

        combined_data = []

        for software, df in software_data.items():
            # 计算准确率
            df = calculate_accuracy(df)

            # 根据基因过滤数据，只选择目标基因
            df_filtered = df[df['Gene'].isin(gene_filter)]

            # 计算每个基因的平均准确率（跨多个样本）
            gene_accuracy = df_filtered.groupby('Gene')['Accuracy'].mean().reset_index()

            # 使用软件名称映射字典，将内部名称替换为期望显示的名称
            gene_accuracy['Software'] = software_name_map.get(software, software)

            # 按照基因分类排序
            gene_accuracy['Gene'] = pd.Categorical(gene_accuracy['Gene'], categories=gene_filter, ordered=True)
            gene_accuracy = gene_accuracy.sort_values('Gene')

            # 将结果添加到 combined_data 列表
            combined_data.append(gene_accuracy)

        # 合并三个软件的数据
        combined_df = pd.concat(combined_data)

        # 创建一个新的图形对象，带有两个子图
        fig, axes = plt.subplots(2, 1, figsize=(fig_width, 3), gridspec_kw={'height_ratios': [2, 1]})

        # ==============================
        # 第一个子图: 准确率柱状图 (Barplot)
        # ==============================
        sns.barplot(x='Gene', y='Accuracy', hue='Software', data=combined_df, palette=colors, ci=None, dodge=True, ax=axes[0])

        axes[0].set_xticklabels([])
        axes[0].set_xlabel('')
        axes[0].set_ylabel('Accuracy', fontsize=12)
        axes[0].set_ylim(0, 1.1)
        axes[0].tick_params(axis='y', labelsize=12)
        axes[0].legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=12, frameon=False)

        # ==============================
        # 第二个子图: 深度分布图 (Raincloud Plot)
        # ==============================
        if depth_data is not None:
            # 如果没有显式的 'Sample' 列，将第一列重命名为 'Sample'
            if 'Sample' not in depth_data.columns:
                depth_data.rename(columns={depth_data.columns[0]: 'Sample'}, inplace=True)

            # 移除基因名称中的 "HLA-" 前缀
            depth_data.columns = depth_data.columns.str.replace('HLA-', '')

            # 将深度数据转换为长格式
            depth_data_melted = depth_data.melt(id_vars=['Sample'], var_name='Gene', value_name='Depth')

            # 根据基因过滤数据，只选择目标基因
            depth_data_melted = depth_data_melted[depth_data_melted['Gene'].isin(gene_filter)]

            # 按照基因分类排序
            depth_data_melted['Gene'] = pd.Categorical(depth_data_melted['Gene'], categories=gene_filter, ordered=True)
            depth_data_melted = depth_data_melted.sort_values('Gene')

            # 设置箱线图的颜色和美化
            boxprops = dict(color='darkgray', linewidth=1, facecolor='#75D5DF', alpha=1)
            medianprops = dict(color='darkgray', linewidth=2)
            whiskerprops = dict(color='darkgray', linewidth=1.5)
            capprops = dict(color='darkgray', linewidth=1.5)

            # 在小提琴图上叠加箱线图（作为雨云图中的“总结统计信息”部分）
            sns.boxplot(
                x='Gene', 
                y='Depth', 
                data=depth_data_melted, 
                width=0.3,  # 控制箱线图宽度
                ax=axes[1],
                showfliers=False,  # 不显示异常值
                boxprops=boxprops,  # 设置箱线图的外观
                medianprops=medianprops,  # 设置中位数线的外观
                whiskerprops=whiskerprops,  # 设置胡须线的外观
                capprops=capprops,  # 设置箱线图顶端和底端线的外观
                zorder=2  # 确保箱线图在小提琴图上方
            )

            # 设置标签和外观
            axes[1].set_xticklabels(gene_filter, rotation=45, ha='right', fontsize=10)
            axes[1].set_xlabel('Locus', fontsize=12)
            axes[1].set_ylabel('Depth', fontsize=12)
            axes[1].tick_params(axis='y', labelsize=12)

        plt.tight_layout()

        # 保存图像为SVG格式
        svg_file_path = f'{dataset}_gene_accuracy_depth_comparison_raincloud.svg'
        plt.savefig(svg_file_path, format='svg', bbox_inches='tight', dpi=600)

        # 显示图像
        plt.show()

# 指定每个数据集的目录
dataset_dirs = {
    'hgsvc2_clr': '../HGSVC2/clr',  
    'hgsvc2_hifi': '../HGSVC2/hifi',  
    'hprc_hifi': '../HPRC/hla/hifi',  
    'hprc_ont': '../HPRC/hla/ont',  
    '1kg_ont': '../1KG_ONT/hla/',
}

# 读取文件并加载数据
all_data = load_data_from_files(dataset_dirs)

# 绘制准确率柱状图并保存为SVG
plot_accuracy(all_data)