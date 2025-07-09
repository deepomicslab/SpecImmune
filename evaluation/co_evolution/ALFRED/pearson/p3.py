#!/usr/bin/env python3
"""
使用dbSNP频率数据计算两个SNP位点ref base ratio的皮尔森相关系数
"""

import numpy as np
from scipy.stats import pearsonr
import pandas as pd

def calculate_correlation_from_dbsnp_data():
    """
    使用dbSNP数据计算相关系数
    """
    
    # 从1000Genomes_30X数据提取（因为这个数据集包含了完整的人群信息）
    # chr1:1379664 (rs2368566) - G>A, ref=G
    rs2368566_data = {
        'Global': 0.2753,      # G频率
        'African': 0.6417,
        'Europe': 0.0869, 
        'South Asian': 0.1772,
        'East Asian': 0.1197,
        'American': 0.157
    }
    
    # chr10:541528 (rs59771674) - A>G, ref=A  
    rs59771674_data = {
        'Global': 0.6918,      # A频率
        'African': 0.3438,
        'Europe': 0.8594,
        'South Asian': 0.8719, 
        'East Asian': 0.7333,
        'American': 0.839
    }
    
    print("=" * 80)
    print("SNP位点信息确认")
    print("=" * 80)
    print(f"位点1: chr1:1379664 → rs2368566 (G>A)")
    print(f"位点2: chr10:541528 → rs59771674 (A>G)")
    print("数据来源: 1000Genomes_30X项目")
    
    print(f"\n{'='*80}")
    print("各人群ref allele频率对比")
    print(f"{'='*80}")
    
    # 使用5个主要人群进行分析（排除Global）
    populations = ['African', 'Europe', 'South Asian', 'East Asian', 'American']
    
    print(f"{'Population':<12} {'rs2368566(G)':<15} {'rs59771674(A)':<15} {'差值':<10}")
    print("-" * 60)
    
    rs2368566_values = []
    rs59771674_values = []
    
    for pop in populations:
        freq1 = rs2368566_data[pop]
        freq2 = rs59771674_data[pop] 
        diff = freq1 - freq2
        
        rs2368566_values.append(freq1)
        rs59771674_values.append(freq2)
        
        print(f"{pop:<12} {freq1:<15.4f} {freq2:<15.4f} {diff:>10.4f}")
    
    # 计算皮尔森相关系数
    correlation, p_value = pearsonr(rs2368566_values, rs59771674_values)
    
    print(f"\n{'='*80}")
    print("皮尔森相关系数分析结果")
    print(f"{'='*80}")
    print(f"相关系数 (r): {correlation:.6f}")
    print(f"P值: {p_value:.6f}")
    print(f"样本数: {len(rs2368566_values)} (人群数)")
    
    # 判断相关性强度和方向
    if abs(correlation) >= 0.8:
        strength = "强"
    elif abs(correlation) >= 0.5:
        strength = "中等"
    elif abs(correlation) >= 0.3:
        strength = "弱"
    else:
        strength = "很弱"
    
    direction = "正" if correlation > 0 else "负"
    significance = "显著" if p_value < 0.05 else "不显著"
    
    print(f"相关性强度: {strength}{direction}相关")
    print(f"统计显著性: {significance} (α=0.05)")
    
    # 详细解释
    print(f"\n{'='*80}")
    print("结果解释")
    print(f"{'='*80}")
    
    if correlation > 0:
        print("正相关表明：当rs2368566的ref allele (G) 频率高的人群中，")
        print("rs59771674的ref allele (A) 频率也倾向于较高")
    else:
        print("负相关表明：当rs2368566的ref allele (G) 频率高的人群中，")
        print("rs59771674的ref allele (A) 频率倾向于较低")
    
    print(f"\nr² = {correlation**2:.4f}，表明两位点ref频率变异的{correlation**2*100:.1f}%")
    print("可以通过线性关系来解释")
    
    # 使用全部数据集验证
    print(f"\n{'='*80}")
    print("其他数据集验证")
    print(f"{'='*80}")
    
    # TopMed数据验证
    topmed_rs2368566 = 0.232207  # G频率
    topmed_rs59771674 = 0.717179  # A频率
    print(f"TopMed数据集:")
    print(f"  rs2368566 ref频率: {topmed_rs2368566:.6f}")
    print(f"  rs59771674 ref频率: {topmed_rs59771674:.6f}")
    print(f"  差值: {topmed_rs2368566 - topmed_rs59771674:.6f}")
    
    return correlation, p_value, rs2368566_values, rs59771674_values

def perform_permutation_test(x, y, n_permutations=10000):
    """
    置换检验验证
    """
    print(f"\n{'='*80}")
    print("置换检验验证")
    print(f"{'='*80}")
    
    observed_corr, _ = pearsonr(x, y)
    perm_corrs = []
    
    np.random.seed(42)
    for i in range(n_permutations):
        y_perm = np.random.permutation(y)
        perm_corr, _ = pearsonr(x, y_perm)
        perm_corrs.append(perm_corr)
    
    perm_corrs = np.array(perm_corrs)
    p_value_perm = np.sum(np.abs(perm_corrs) >= np.abs(observed_corr)) / n_permutations
    
    print(f"置换检验P值: {p_value_perm:.6f}")
    print(f"95%置信区间: [{np.percentile(perm_corrs, 2.5):.4f}, {np.percentile(perm_corrs, 97.5):.4f}]")
    
    return p_value_perm

def main():
    """
    主函数
    """
    try:
        correlation, p_value, x, y = calculate_correlation_from_dbsnp_data()
        p_value_perm = perform_permutation_test(x, y)
        
        print(f"\n{'='*80}")
        print("最终结论")
        print(f"{'='*80}")
        print(f"rs2368566 (chr1:1379664) 和 rs59771674 (chr10:541528)")
        print(f"两位点ref allele频率的皮尔森相关系数为: r = {correlation:.4f}")
        
        if p_value < 0.05 and p_value_perm < 0.05:
            print("相关性在统计学上显著 (p < 0.05)")
        else:
            print("相关性在统计学上不显著 (p ≥ 0.05)")
            
        return 0
        
    except Exception as e:
        print(f"计算过程中出现错误: {e}")
        return 1

if __name__ == "__main__":
    exit(main())