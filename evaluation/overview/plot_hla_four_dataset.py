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
genes_for_hlala = ['A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRA', 'DRB1', 'E', 'F', 'G', 'H', 'K', 'V', 'DRB3', 'DRB4']
gene_order = genes_for_hlala

KGP_genes = ['A', 'B', 'C', 'DQB1', 'DRB1']
colors = ["#D2D2D2", "#062565", "#0098B4"]
colors = ["#98B6C4", "#8D93AF", "#3E4271"]

outdir="revise_as_r1"
if not os.path.exists(outdir):
    os.makedirs(outdir)


# 软件名称映射字典
software_name_map = {
    'hlala': 'HLA*LA',
    'spechla': 'SpecHLA',
    'speclong': 'SpecImmune'
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
            if dataset == '1kg_ont' and software == 'hlala':
                file_path = os.path.join(data_dir, f'{software}.match.3.csv')
            if dataset == 'hgsvc2_clr' and software == 'hlala':
                file_path = os.path.join(data_dir, f'{software}.match.qv.csv')
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
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# 绘图函数，只包含准确率柱状图
def plot_accuracy(all_data):
    for dataset, data in all_data.items():
        software_data = data['software_data']

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

        # 创建单个图形对象
        fig, ax = plt.subplots(1, 1, figsize=(fig_width, 2.5))

        # ==============================
        # 准确率柱状图 (Barplot)
        # ==============================
        sns.barplot(x='Gene', y='Accuracy', hue='Software', data=combined_df, palette=colors, ci=None, dodge=True, ax=ax)

        ax.set_xticklabels(gene_filter, rotation=45, ha='right', fontsize=10)
        ax.set_xlabel('Locus', fontsize=12)
        ax.set_ylabel('Accuracy', fontsize=12)
        ax.set_ylim(0, 1.1)
        ax.tick_params(axis='y', labelsize=12)
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=12, frameon=False)
        
        # ========== 添加淡色横向网格线 ==========
        ax.yaxis.grid(True, linestyle='--', alpha=0.3, linewidth=0.8, color='gray')
        ax.set_axisbelow(True)  # 将网格线放在图形元素下方
        # =======================================

        plt.tight_layout()

        # 保存图像为SVG格式
        svg_file_path = f'{outdir}/{dataset}_gene_accuracy_comparison.svg'
        plt.savefig(svg_file_path, format='svg', bbox_inches='tight', dpi=600)

        # 关闭图形以释放内存
        plt.close()
        # 显示图像
        # plt.show()

# 指定每个数据集的目录
dataset_dirs = {
    'hgsvc2_clr': '../HGSVC2/clr',  
    'hgsvc2_hifi': '../HGSVC2/hifi',  
    'hprc_hifi': '../HPRC/hla/hifi',  
    'hprc_ont': '../HPRC/hla/ont',  
    '1kg_ont': '../1KG_ONT/hla/',
}
# dataset_dirs = { 
#     'hgsvc2_clr': '../HGSVC2/clr', 
# }

# 读取文件并加载数据
all_data = load_data_from_files(dataset_dirs)

# 绘制准确率柱状图并保存为SVG
plot_accuracy(all_data)