import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# 定义数据集目录
dataset_dirs = {
    'hgsvc2_clr': '../HGSVC2/clr',  
    'hgsvc2_hifi': '../HGSVC2/hifi',  
    'hprc_hifi': '../HPRC/hla/hifi',  
    'hprc_ont': '../HPRC/hla/ont',  
}

# 定义深度阈值
depth_cutoffs = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50]

# 定义每个软件需要考虑的基因
genes_for_hlala = ['A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRA', 'DRB1', 'E', 'F', 'G', 'H', 'K', 'V', 'DRB3', 'DRB4']
genes_for_spechla = ['A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRB1']

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

# 计算基因类别的准确率
def calculate_accuracies_by_class(df, gene_reads_df, software_genes, chosen_depth):
    class_level_data = {class_name: {'<': [], '>=': []} for class_name in gene_classes}
    match_totals = {class_name: {'<': {'Match': 0, 'Total': 0}, '>=': {'Match': 0, 'Total': 0}} for class_name in gene_classes}

    relevant_genes = software_genes
    df = df[df['Gene'].isin(relevant_genes)]

    for gene in relevant_genes:
        samples_lt_chosen = gene_reads_df[gene][gene_reads_df[gene] < chosen_depth].dropna().index
        samples_ge_chosen = gene_reads_df[gene][gene_reads_df[gene] >= chosen_depth].dropna().index

        filtered_df_lt_chosen = df[(df['Gene'] == gene) & (df['Sample'].isin(samples_lt_chosen))]
        filtered_df_ge_chosen = df[(df['Gene'] == gene) & (df['Sample'].isin(samples_ge_chosen))]

        # 计算 <chosen_depth 的准确率
        if not filtered_df_lt_chosen.empty:
            total_match_lt_chosen = filtered_df_lt_chosen['Match'].sum()
            total_total_lt_chosen = filtered_df_lt_chosen['Total'].sum()
            accuracy_lt_chosen = total_match_lt_chosen / total_total_lt_chosen if total_total_lt_chosen > 0 else None
        else:
            total_match_lt_chosen = 0
            total_total_lt_chosen = 0
            accuracy_lt_chosen = None

        # 计算 >=chosen_depth 的准确率
        if not filtered_df_ge_chosen.empty:
            total_match_ge_chosen = filtered_df_ge_chosen['Match'].sum()
            total_total_ge_chosen = filtered_df_ge_chosen['Total'].sum()
            accuracy_ge_chosen = total_match_ge_chosen / total_total_ge_chosen if total_total_ge_chosen > 0 else None
        else:
            total_match_ge_chosen = 0
            total_total_ge_chosen = 0
            accuracy_ge_chosen = None

        # 按基因类别保存准确率和匹配/总数
        for class_name, class_genes in gene_classes.items():
            if gene in class_genes:
                class_level_data[class_name]['<'].append(accuracy_lt_chosen)
                class_level_data[class_name]['>='].append(accuracy_ge_chosen)
                match_totals[class_name]['<']['Match'] += total_match_lt_chosen
                match_totals[class_name]['<']['Total'] += total_total_lt_chosen
                match_totals[class_name]['>=']['Match'] += total_match_ge_chosen
                match_totals[class_name]['>=']['Total'] += total_total_ge_chosen

    return class_level_data, match_totals

# 用于保存输出到 CSV 的数据
csv_data = []

# 遍历每个数据集并计算每个深度阈值下的准确率
for dataset_name, dataset_dir in dataset_dirs.items():
    print(f"\nDataset: {dataset_name}")

    # 加载这个数据集的 CSV 文件
    df1 = pd.read_csv(f'{dataset_dir}/hlala.match.csv')
    df2 = pd.read_csv(f'{dataset_dir}/spechla.match.csv')
    df3 = pd.read_csv(f'{dataset_dir}/speclong.match.csv')

    # 定义 genes_for_speclong
    genes_for_speclong = [gene for gene in pd.concat([df1['Gene'], df2['Gene'], df3['Gene']]).unique() if gene not in ['']]

    # 加载基因读深度
    gene_reads_df = pd.read_csv(f'{dataset_dir}/read_depth.csv', index_col=0)
    gene_reads_df.columns = gene_reads_df.columns.str.replace('HLA-', '', regex=False)

    dfs = [df1, df2, df3]
    df_names = ['HLA*LA', 'SpecHLA', 'SpecLong']

    software_genes = {
        'HLA*LA': genes_for_hlala,
        'SpecHLA': genes_for_spechla,
        'SpecLong': genes_for_speclong
    }

    # 遍历每个深度阈值
    for chosen_depth in depth_cutoffs:
        print(f"\nDepth cutoff: {chosen_depth}")

        for df, df_name in zip(dfs, df_names):
            # 计算每个基因类别的准确率和匹配/总数
            class_level_data, match_totals = calculate_accuracies_by_class(df, gene_reads_df, software_genes[df_name], chosen_depth)

            # 计算并收集每个基因类别的平均准确率和匹配/总数
            for class_name, accuracies in class_level_data.items():
                accuracies_lt_chosen = [acc for acc in accuracies['<'] if acc is not None]
                accuracies_ge_chosen = [acc for acc in accuracies['>='] if acc is not None]

                average_accuracy_lt_chosen = sum(accuracies_lt_chosen) / len(accuracies_lt_chosen) if accuracies_lt_chosen else None
                average_accuracy_ge_chosen = sum(accuracies_ge_chosen) / len(accuracies_ge_chosen) if accuracies_ge_chosen else None

                # 获取匹配/总数
                total_match_lt = match_totals[class_name]['<']['Match']
                total_total_lt = match_totals[class_name]['<']['Total']
                total_match_ge = match_totals[class_name]['>=']['Match']
                total_total_ge = match_totals[class_name]['>=']['Total']

                # 保存数据到列表
                csv_data.append({
                    'Dataset': dataset_name,
                    'Software': df_name,
                    'Class': class_name,
                    'Depth Cutoff': chosen_depth,
                    'Average Accuracy < Depth': average_accuracy_lt_chosen,
                    'Match < Depth': total_match_lt,
                    'Total < Depth': total_total_lt,
                    'Average Accuracy >= Depth': average_accuracy_ge_chosen,
                    'Match >= Depth': total_match_ge,
                    'Total >= Depth': total_total_ge
                })

# 将数据转换为 DataFrame
csv_df = pd.DataFrame(csv_data)

# 将 DataFrame 保存为 CSV 文件
csv_df.to_csv('average_accuracies_by_class_per_depth.csv', index=False)

print("CSV 文件 'average_accuracies_by_class_per_depth.csv' 已创建。")