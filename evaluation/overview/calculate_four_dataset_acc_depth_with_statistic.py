import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests
import warnings
warnings.filterwarnings('ignore')

# ==================== 配置部分 ====================
# 定义数据集目录
dataset_dirs = {
    'hgsvc2_clr': '../HGSVC2/clr',  
    'hgsvc2_hifi': '../HGSVC2/hifi',  
    'hprc_hifi': '../HPRC/hla/hifi',  
    'hprc_ont': '../HPRC/hla/ont',
    '1kg_ont': '../1KG_ONT/hla',  # 添加 1kg_ont 数据集
}


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

# ==================== 函数定义 ====================

def classify_gene(gene_name):
    """判断基因属于哪个类别"""
    for class_name, genes in gene_classes.items():
        if gene_name in genes:
            return class_name
    return 'Other'

def calculate_typing_success_rate(df_list, df_names, genes_dict):
    """
    计算每个软件的分型成功率
    """
    results = []
    
    for df, software_name in zip(df_list, df_names):
        # 筛选该软件关注的基因
        relevant_genes = genes_dict[software_name]
        df_filtered = df[df['Gene'].isin(relevant_genes)].copy()
        
        # 添加基因分类
        df_filtered['Class'] = df_filtered['Gene'].apply(classify_gene)
        
        # 按类别统计
        for class_name in list(gene_classes.keys()) + ['All Genes']:
            if class_name == 'All Genes':
                class_data = df_filtered
            else:
                class_data = df_filtered[df_filtered['Class'] == class_name]
            
            if len(class_data) == 0:
                continue
            
            # 总样本数（基因-样本对）
            total_samples = len(class_data)
            
            # 成功分型的样本数（Total > 0 表示有分型结果）
            typed_samples = len(class_data[class_data['Total'] > 0])
            
            typing_rate = typed_samples / total_samples if total_samples > 0 else 0
            
            results.append({
                'Software': software_name,
                'Class': class_name,
                'Total_Samples': total_samples,
                'Typed_Samples': typed_samples,
                'Typing_Success_Rate': typing_rate
            })
    
    return pd.DataFrame(results)

def calculate_accuracy_common_samples(df_list, df_names, genes_dict):
    """
    计算在所有软件共同成功分型的样本上的准确率
    """
    results = []
    
    # 首先，需要合并所有软件的数据以找到共同样本
    df_merged = None
    
    for idx, (df, software_name) in enumerate(zip(df_list, df_names)):
        relevant_genes = genes_dict[software_name]
        df_filtered = df[df['Gene'].isin(relevant_genes)].copy()
        df_filtered['Class'] = df_filtered['Gene'].apply(classify_gene)
        
        # 创建唯一标识
        df_filtered['ID'] = df_filtered['Gene'].astype(str) + '_' + df_filtered['Sample'].astype(str)
        
        # 重命名列
        df_filtered = df_filtered.rename(columns={
            'Match': f'{software_name}_Match',
            'Total': f'{software_name}_Total'
        })
        
        if df_merged is None:
            df_merged = df_filtered[['ID', 'Gene', 'Sample', 'Class', 
                                      f'{software_name}_Match', f'{software_name}_Total']].copy()
        else:
            df_merged = df_merged.merge(
                df_filtered[['ID', f'{software_name}_Match', f'{software_name}_Total']],
                on='ID',
                how='outer'
            )
    
    # 填充缺失值为0（表示该软件没有对这个基因-样本对进行分型）
    for software_name in df_names:
        df_merged[f'{software_name}_Match'] = df_merged[f'{software_name}_Match'].fillna(0)
        df_merged[f'{software_name}_Total'] = df_merged[f'{software_name}_Total'].fillna(0)
    
    # 计算每个软件在共同样本上的准确率
    for software_name in df_names:
        # 按类别分组
        for class_name in list(gene_classes.keys()) + ['All Genes']:
            if class_name == 'All Genes':
                class_data = df_merged.copy()
            else:
                class_data = df_merged[df_merged['Class'] == class_name].copy()
            
            if len(class_data) == 0:
                continue
            
            # 筛选所有软件都有分型结果的样本
            # 使用布尔索引，确保索引对齐
            common_mask = pd.Series([True] * len(class_data), index=class_data.index)
            for sw in df_names:
                common_mask = common_mask & (class_data[f'{sw}_Total'] > 0)
            
            common_data = class_data[common_mask].copy()
            
            if len(common_data) == 0:
                continue
            
            # 计算该软件的准确率
            total_match = common_data[f'{software_name}_Match'].sum()
            total_alleles = common_data[f'{software_name}_Total'].sum()
            
            accuracy = total_match / total_alleles if total_alleles > 0 else 0
            
            # 计算置信区间 (Wilson score)
            if total_alleles > 0:
                p = accuracy
                n = total_alleles
                z = 1.96
                denominator = 1 + z**2/n
                centre = (p + z**2/(2*n)) / denominator
                interval = z * np.sqrt(p*(1-p)/n + z**2/(4*n**2)) / denominator
                ci_low = max(0, centre - interval)
                ci_upp = min(1, centre + interval)
            else:
                ci_low, ci_upp = 0, 0
            
            results.append({
                'Software': software_name,
                'Class': class_name,
                'Common_Samples': len(common_data),
                'Total_Matches': int(total_match),
                'Total_Alleles': int(total_alleles),
                'Accuracy_Common': accuracy,
                'CI_Low': ci_low,
                'CI_Upp': ci_upp
            })
    
    return pd.DataFrame(results), df_merged

def calculate_accuracy_with_penalty(df_list, df_names, genes_dict):
    """
    计算包含分型失败惩罚的准确率
    未分型样本（Total=0）计为错误
    """
    results = []
    
    for df, software_name in zip(df_list, df_names):
        relevant_genes = genes_dict[software_name]
        df_filtered = df[df['Gene'].isin(relevant_genes)].copy()
        df_filtered['Class'] = df_filtered['Gene'].apply(classify_gene)
        
        for class_name in list(gene_classes.keys()) + ['All Genes']:
            if class_name == 'All Genes':
                class_data = df_filtered
            else:
                class_data = df_filtered[df_filtered['Class'] == class_name]
            
            if len(class_data) == 0:
                continue
            
            # 这里假设每个基因-样本对应该有2个等位基因
            # 如果实际情况不同，需要调整
            expected_total = len(class_data) * 2  # 假设每个样本2个等位基因
            
            # 实际匹配数
            actual_match = class_data['Match'].sum()
            
            # 实际分型的等位基因总数
            actual_total = class_data['Total'].sum()
            
            # 包含惩罚的准确率
            accuracy_penalty = actual_match / expected_total if expected_total > 0 else 0
            
            # 不含惩罚的准确率（仅成功分型的）
            accuracy_no_penalty = actual_match / actual_total if actual_total > 0 else 0
            
            results.append({
                'Software': software_name,
                'Class': class_name,
                'Expected_Total_Alleles': expected_total,
                'Actual_Total_Alleles': int(actual_total),
                'Total_Matches': int(actual_match),
                'Accuracy_With_Penalty': accuracy_penalty,
                'Accuracy_No_Penalty': accuracy_no_penalty
            })
    
    return pd.DataFrame(results)

def perform_statistical_tests(df_merged, df_names, software_pairs):
    """
    执行统计检验
    """
    results = []
    
    for sw1, sw2 in software_pairs:
        for class_name in list(gene_classes.keys()) + ['All Genes']:
            if class_name == 'All Genes':
                class_data = df_merged.copy()
            else:
                class_data = df_merged[df_merged['Class'] == class_name].copy()
            
            if len(class_data) == 0:
                continue
            
            # 筛选两个软件都有分型结果的样本
            valid_mask = (class_data[f'{sw1}_Total'] > 0) & (class_data[f'{sw2}_Total'] > 0)
            valid_data = class_data[valid_mask].copy()
            
            if len(valid_data) == 0:
                continue
            
            # 计算准确率
            sw1_match = valid_data[f'{sw1}_Match'].sum()
            sw1_total = valid_data[f'{sw1}_Total'].sum()
            sw2_match = valid_data[f'{sw2}_Match'].sum()
            sw2_total = valid_data[f'{sw2}_Total'].sum()
            
            sw1_acc = sw1_match / sw1_total if sw1_total > 0 else 0
            sw2_acc = sw2_match / sw2_total if sw2_total > 0 else 0
            
            # Fisher精确检验
            contingency = [
                [int(sw1_match), int(sw1_total - sw1_match)],
                [int(sw2_match), int(sw2_total - sw2_match)]
            ]
            
            try:
                fisher_p = stats.fisher_exact(contingency)[1]
            except:
                fisher_p = 1.0
            
            # McNemar检验 (样本级别)
            both_correct = 0
            sw1_only_correct = 0
            sw2_only_correct = 0
            both_wrong = 0
            
            for _, row in valid_data.iterrows():
                sw1_correct = (row[f'{sw1}_Match'] == row[f'{sw1}_Total']) and (row[f'{sw1}_Total'] > 0)
                sw2_correct = (row[f'{sw2}_Match'] == row[f'{sw2}_Total']) and (row[f'{sw2}_Total'] > 0)
                
                if sw1_correct and sw2_correct:
                    both_correct += 1
                elif sw1_correct and not sw2_correct:
                    sw1_only_correct += 1
                elif not sw1_correct and sw2_correct:
                    sw2_only_correct += 1
                else:
                    both_wrong += 1
            
            # McNemar检验
            try:
                if sw1_only_correct + sw2_only_correct > 0:
                    mcnemar_p = stats.binom_test(
                        min(sw1_only_correct, sw2_only_correct),
                        sw1_only_correct + sw2_only_correct,
                        p=0.5
                    )
                else:
                    mcnemar_p = 1.0
            except:
                mcnemar_p = 1.0
            
            # 效应量
            diff = abs(sw1_acc - sw2_acc)
            if diff < 0.02:
                effect = 'Negligible'
            elif diff < 0.05:
                effect = 'Small'
            elif diff < 0.10:
                effect = 'Medium'
            else:
                effect = 'Large'
            
            results.append({
                'Comparison': f'{sw1} vs {sw2}',
                'Class': class_name,
                'Common_Samples': len(valid_data),
                f'{sw1}_Accuracy': sw1_acc,
                f'{sw2}_Accuracy': sw2_acc,
                'Accuracy_Difference': sw1_acc - sw2_acc,
                'Fisher_P': fisher_p,
                'McNemar_P': mcnemar_p,
                'Effect_Size': effect,
                'Both_Correct': both_correct,
                f'{sw1}_Only_Correct': sw1_only_correct,
                f'{sw2}_Only_Correct': sw2_only_correct,
                'Both_Wrong': both_wrong
            })
    
    df_results = pd.DataFrame(results)
    
    # 多重检验校正
    if len(df_results) > 0:
        fisher_corrected = multipletests(df_results['Fisher_P'], method='fdr_bh')
        mcnemar_corrected = multipletests(df_results['McNemar_P'], method='fdr_bh')
        
        df_results['Fisher_P_Adjusted'] = fisher_corrected[1]
        df_results['Fisher_Significant'] = fisher_corrected[0]
        df_results['McNemar_P_Adjusted'] = mcnemar_corrected[1]
        df_results['McNemar_Significant'] = mcnemar_corrected[0]
    
    return df_results

# ==================== 主程序 ====================

# 存储所有数据集的结果
all_typing_rates = []
all_accuracy_common = []
all_accuracy_penalty = []
all_statistics = []

# 遍历每个数据集
for dataset_name, dataset_dir in dataset_dirs.items():
    print(f"\n{'='*80}")
    print(f"处理数据集: {dataset_name}")
    print(f"{'='*80}")

    # 加载数据
    try:
        # 根据数据集名称选择正确的文件名
        if dataset_name == 'hgsvc2_clr':
            hlala_file = f'{dataset_dir}/hlala.match.qv.csv'
        elif dataset_name == '1kg_ont':
            hlala_file = f'{dataset_dir}/hlala.match.3.csv'
        else:
            hlala_file = f'{dataset_dir}/hlala.match.csv'
        
        df1 = pd.read_csv(hlala_file)
        df2 = pd.read_csv(f'{dataset_dir}/spechla.match.csv')
        df3 = pd.read_csv(f'{dataset_dir}/speclong.match.csv')
        
        print(f"✓ 成功加载数据文件")
        print(f"  - HLA*LA: {hlala_file}")
        print(f"  - SpecHLA: {dataset_dir}/spechla.match.csv")
        print(f"  - SpecLong: {dataset_dir}/speclong.match.csv")
        
    except FileNotFoundError as e:
        print(f"⚠ 警告: 文件未找到 - {e}")
        print(f"  跳过数据集: {dataset_name}")
        continue

    # 定义 genes_for_speclong
    genes_for_speclong = [gene for gene in pd.concat([df1['Gene'], df2['Gene'], df3['Gene']]).unique() if gene not in ['', None]]

    dfs = [df1, df2, df3]
    df_names = ['HLA*LA', 'SpecHLA', 'SpecLong']

    software_genes = {
        'HLA*LA': genes_for_hlala,
        'SpecHLA': genes_for_spechla,
        'SpecLong': genes_for_speclong
    }

    # 1. 计算分型成功率
    print("\n【一、分型成功率】")
    df_typing = calculate_typing_success_rate(dfs, df_names, software_genes)
    df_typing['Dataset'] = dataset_name
    all_typing_rates.append(df_typing)
    
    for software in df_names:
        sw_data = df_typing[df_typing['Software'] == software]
        for class_name in ['Class I Genes', 'Class II Genes', 'All Genes']:
            class_data = sw_data[sw_data['Class'] == class_name]
            if len(class_data) > 0:
                row = class_data.iloc[0]
                print(f"{software} - {class_name}: {row['Typing_Success_Rate']*100:.2f}% "
                      f"({row['Typed_Samples']}/{row['Total_Samples']})")

    # 2. 计算共同样本准确率
    print("\n【二、准确率（共同样本）】")
    df_acc_common, df_merged = calculate_accuracy_common_samples(dfs, df_names, software_genes)
    df_acc_common['Dataset'] = dataset_name
    all_accuracy_common.append(df_acc_common)
    
    for class_name in ['Class I Genes', 'Class II Genes', 'All Genes']:
        print(f"\n{class_name}:")
        class_data = df_acc_common[df_acc_common['Class'] == class_name]
        for _, row in class_data.iterrows():
            print(f"  {row['Software']}: {row['Accuracy_Common']*100:.2f}% "
                  f"[{row['CI_Low']*100:.2f}%-{row['CI_Upp']*100:.2f}%] "
                  f"({row['Total_Matches']}/{row['Total_Alleles']})")

    # 3. 计算含惩罚准确率
    print("\n【三、准确率（含分型失败惩罚）】")
    df_acc_penalty = calculate_accuracy_with_penalty(dfs, df_names, software_genes)
    df_acc_penalty['Dataset'] = dataset_name
    all_accuracy_penalty.append(df_acc_penalty)
    
    for class_name in ['Class I Genes', 'Class II Genes', 'All Genes']:
        print(f"\n{class_name}:")
        class_data = df_acc_penalty[df_acc_penalty['Class'] == class_name]
        for _, row in class_data.iterrows():
            print(f"  {row['Software']}: {row['Accuracy_With_Penalty']*100:.2f}% "
                  f"({row['Total_Matches']}/{row['Expected_Total_Alleles']})")

    # 4. 统计检验
    print("\n【四、统计显著性检验】")
    software_pairs = [
        ('SpecLong', 'HLA*LA'),
        ('SpecLong', 'SpecHLA'),
        ('HLA*LA', 'SpecHLA')
    ]
    
    df_stats = perform_statistical_tests(df_merged, df_names, software_pairs)
    df_stats['Dataset'] = dataset_name
    all_statistics.append(df_stats)
    
    significant = df_stats[df_stats['Fisher_Significant'] == True]
    print(f"\n发现 {len(significant)} 项显著差异:")
    
    for _, row in significant.iterrows():
        print(f"\n{row['Comparison']} - {row['Class']}:")
        print(f"  准确率差异: {row['Accuracy_Difference']*100:+.2f}%")
        print(f"  效应量: {row['Effect_Size']}")
        print(f"  Fisher p (校正): {row['Fisher_P_Adjusted']:.4f}")
        print(f"  McNemar p (校正): {row['McNemar_P_Adjusted']:.4f}")

# 合并所有数据集的结果
print(f"\n\n{'='*80}")
print("保存结果文件...")
print(f"{'='*80}")

df_all_typing = pd.concat(all_typing_rates, ignore_index=True)
df_all_acc_common = pd.concat(all_accuracy_common, ignore_index=True)
df_all_acc_penalty = pd.concat(all_accuracy_penalty, ignore_index=True)
df_all_stats = pd.concat(all_statistics, ignore_index=True)

# 保存CSV文件
output_files = {
    'typing_success_rates.csv': df_all_typing,
    'accuracy_common_samples.csv': df_all_acc_common,
    'accuracy_with_penalty.csv': df_all_acc_penalty,
    'statistical_comparisons.csv': df_all_stats
}

for filename, df in output_files.items():
    df.to_csv(filename, index=False)
    print(f"✓ 已保存: {filename}")

print("\n" + "="*80)
print("分析完成！")
print("="*80)
print(f"\n处理的数据集: {', '.join(dataset_dirs.keys())}")
print(f"生成的输出文件: {len(output_files)} 个")