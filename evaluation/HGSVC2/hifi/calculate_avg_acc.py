import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Step 1: 读取 CSV 文件到 DataFrame 中
df1 = pd.read_csv('hlala.match.csv')
df2 = pd.read_csv('spechla.match.csv')
df3 = pd.read_csv('speclong.match.csv')

# 将所有 DataFrame 放入列表中，同时保存对应的软件名称
dfs = [df1, df2, df3]
df_names = ['HLA*LA', 'SpecHLA', 'SpecLong']  # 软件名称

# Step 2: 定义各软件需要考虑的基因列表
genes_for_hlala = [
    'A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1',
    'DQB1', 'DRA', 'DRB1', 'E', 'F', 'G',
    'H', 'K', 'V', 'DRB3', 'DRB4'
]
genes_for_spechla = ['A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRB1']
genes_for_speclong = [
    gene for gene in pd.concat([df1['Gene'], df2['Gene'], df3['Gene']]).unique()
    if gene not in ['HFE', 'DRB3', 'MICA', 'DPA2', 'MICB', 'DOA', 'W', 'U', 'DOA2', 'S']
]

# Step 3: 读取基因读深数据，并清洗基因名称（去掉 "HLA-" 前缀）
gene_reads_df = pd.read_csv('read_depth.csv', index_col=0)
gene_reads_df.columns = gene_reads_df.columns.str.replace('HLA-', '', regex=False)

# 建立字典，将软件名称映射到各自对应的基因列表
software_genes = {
    'HLA*LA': genes_for_hlala,
    'SpecHLA': genes_for_spechla,
    'SpecLong': genes_for_speclong
}

# Step 4: 设置阈值
chosen_depth = 10  # 例如：10

# Step 5: 初始化字典，用于存储各深度类别的处理结果
# 每个记录中存储：
#   'Gene'         : 基因名称
#   'Accuracy'     : 经过过滤后的准确率（只针对有效数据）
#   'Total_Match'  : 经过过滤后的正确等位基因数
#   'Total_Total'  : 经过过滤后的总等位基因数
plot_data = {
    f'<{chosen_depth}': {df_name: [] for df_name in df_names},
    f'>={chosen_depth}': {df_name: [] for df_name in df_names}
}

# Step 6: 对每个 DataFrame 中的每个基因计算经过阈值过滤后的准确率及计数
for df, df_name in zip(dfs, df_names):
    # 只对该软件关注的基因进行处理
    relevant_genes = software_genes[df_name]
    df = df[df['Gene'].isin(relevant_genes)]
    
    for gene in relevant_genes:
        # 从 gene_reads_df 中获取该基因读深对应的样本索引
        samples_lt = gene_reads_df[gene][gene_reads_df[gene] < chosen_depth].dropna().index
        samples_ge = gene_reads_df[gene][gene_reads_df[gene] >= chosen_depth].dropna().index
        
        # 处理读深 < chosen_depth 的情况
        filtered_df_lt = df[(df['Gene'] == gene) & (df['Sample'].isin(samples_lt))]
        if not filtered_df_lt.empty:
            total_match_lt = filtered_df_lt['Match'].sum()
            total_total_lt = filtered_df_lt['Total'].sum()
            accuracy_lt = total_match_lt / total_total_lt if total_total_lt > 0 else None
        else:
            accuracy_lt = None
            total_match_lt = 0
            total_total_lt = 0
        
        plot_data[f'<{chosen_depth}'][df_name].append({
            'Gene': gene,
            'Accuracy': accuracy_lt,
            'Total_Match': total_match_lt,
            'Total_Total': total_total_lt
        })
        
        # 处理读深 >= chosen_depth 的情况
        filtered_df_ge = df[(df['Gene'] == gene) & (df['Sample'].isin(samples_ge))]
        if not filtered_df_ge.empty:
            total_match_ge = filtered_df_ge['Match'].sum()
            total_total_ge = filtered_df_ge['Total'].sum()
            accuracy_ge = total_match_ge / total_total_ge if total_total_ge > 0 else None
        else:
            accuracy_ge = None
            total_match_ge = 0
            total_total_ge = 0
        
        plot_data[f'>={chosen_depth}'][df_name].append({
            'Gene': gene,
            'Accuracy': accuracy_ge,
            'Total_Match': total_match_ge,
            'Total_Total': total_total_ge
        })

# Step 7: 输出结果，只统计经过阈值过滤且有有效数据（Accuracy 非 None）的基因
for depth_category in plot_data:
    print(f"\nResults for depth category: {depth_category}")
    for df_name in df_names:
        print(f"\nSoftware: {df_name}")
        
        # 仅使用 Accuracy 有效的条目
        valid_entries = [
            entry for entry in plot_data[depth_category][df_name]
            if entry['Accuracy'] is not None
        ]
        
        if not valid_entries:
            print("  No valid gene data for evaluation")
            continue
        
        accuracies = []
        total_matches_sum = 0  # 仅累计有效基因的正确等位基因数
        total_alleles_sum = 0  # 仅累计有效基因的总等位基因数
        
        for entry in valid_entries:
            gene = entry['Gene']
            accuracy = entry['Accuracy']
            accuracies.append(accuracy)
            total_matches_sum += entry['Total_Match']
            total_alleles_sum += entry['Total_Total']
            
            if entry['Total_Total'] > 0:
                print(f"  Gene: {gene}, Accuracy: {accuracy:.2f}, Correct Alleles: {entry['Total_Match']}, Total Alleles: {entry['Total_Total']}")
            else:
                print(f"  Gene: {gene}, Accuracy: No valid data")
        
        average_accuracy = sum(accuracies) / len(accuracies)
        print(f"\n  Total Correct Alleles: {total_matches_sum}, Total Processed Alleles: {total_alleles_sum}")
        print(f"  Average Accuracy for {df_name}: {average_accuracy:.2f}")