#!/usr/bin/env python3
"""
计算两个SNP位点ref base ratio的皮尔森相关系数
包含置换检验验证统计显著性
"""

import numpy as np
from scipy.stats import pearsonr
import random
import matplotlib.pyplot as plt
import seaborn as sns

def permutation_test(x, y, n_permutations=10000, random_seed=42):
    """
    执行置换检验
    
    参数:
    x, y: 两组数据
    n_permutations: 置换次数
    random_seed: 随机种子
    
    返回:
    observed_corr: 观察到的相关系数
    perm_corrs: 置换相关系数数组
    p_value_perm: 置换检验p值
    """
    # 设置随机种子以确保结果可重现
    random.seed(random_seed)
    np.random.seed(random_seed)
    
    # 计算观察到的相关系数
    observed_corr, _ = pearsonr(x, y)
    
    # 存储置换结果
    perm_corrs = []
    
    print(f"正在进行置换检验（{n_permutations}次置换）...")
    
    for i in range(n_permutations):
        if (i + 1) % 1000 == 0:
            print(f"已完成 {i + 1}/{n_permutations} 次置换")
        
        # 随机打乱其中一组数据
        y_permuted = np.random.permutation(y)
        
        # 计算置换后的相关系数
        perm_corr, _ = pearsonr(x, y_permuted)
        perm_corrs.append(perm_corr)
    
    perm_corrs = np.array(perm_corrs)
    
    # 计算双尾检验的p值
    # 计算极端值的比例（绝对值大于等于观察值绝对值的比例）
    extreme_count = np.sum(np.abs(perm_corrs) >= np.abs(observed_corr))
    p_value_perm = extreme_count / n_permutations
    
    return observed_corr, perm_corrs, p_value_perm

def plot_permutation_results(observed_corr, perm_corrs, p_value_perm, save_plot=False):
    """
    绘制置换检验结果图
    """
    plt.figure(figsize=(12, 8))
    
    # 创建子图
    plt.subplot(2, 2, 1)
    # 绘制置换相关系数的直方图
    plt.hist(perm_corrs, bins=50, alpha=0.7, color='skyblue', edgecolor='black')
    plt.axvline(observed_corr, color='red', linestyle='--', linewidth=2, 
                label=f'观察值 r = {observed_corr:.4f}')
    plt.axvline(-observed_corr, color='red', linestyle='--', linewidth=2, alpha=0.5)
    plt.xlabel('置换相关系数')
    plt.ylabel('频次')
    plt.title(f'置换检验结果分布\nP值 = {p_value_perm:.4f}')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Q-Q图检验正态性
    plt.subplot(2, 2, 2)
    from scipy import stats
    stats.probplot(perm_corrs, dist="norm", plot=plt)
    plt.title('置换结果Q-Q图')
    plt.grid(True, alpha=0.3)
    
    # 累积分布图
    plt.subplot(2, 2, 3)
    sorted_perm = np.sort(perm_corrs)
    cumulative = np.arange(1, len(sorted_perm) + 1) / len(sorted_perm)
    plt.plot(sorted_perm, cumulative, 'b-', linewidth=2)
    plt.axvline(observed_corr, color='red', linestyle='--', linewidth=2,
                label=f'观察值 r = {observed_corr:.4f}')
    plt.xlabel('相关系数')
    plt.ylabel('累积概率')
    plt.title('置换结果累积分布')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # 统计摘要
    plt.subplot(2, 2, 4)
    plt.axis('off')
    
    # 计算统计量
    mean_perm = np.mean(perm_corrs)
    std_perm = np.std(perm_corrs)
    min_perm = np.min(perm_corrs)
    max_perm = np.max(perm_corrs)
    
    # 计算置信区间
    ci_95 = np.percentile(perm_corrs, [2.5, 97.5])
    ci_99 = np.percentile(perm_corrs, [0.5, 99.5])
    
    stats_text = f"""置换检验统计摘要:
    
观察相关系数: {observed_corr:.6f}
置换检验P值: {p_value_perm:.6f}

置换结果统计:
均值: {mean_perm:.6f}
标准差: {std_perm:.6f}
最小值: {min_perm:.6f}
最大值: {max_perm:.6f}

95%置信区间: [{ci_95[0]:.4f}, {ci_95[1]:.4f}]
99%置信区间: [{ci_99[0]:.4f}, {ci_99[1]:.4f}]

极端值数量: {np.sum(np.abs(perm_corrs) >= np.abs(observed_corr))}
总置换次数: {len(perm_corrs)}"""
    
    plt.text(0.05, 0.95, stats_text, transform=plt.gca().transAxes, 
             fontsize=10, verticalalignment='top', fontfamily='monospace')
    
    plt.tight_layout()
    
    if save_plot:
        plt.savefig('permutation_test_results.png', dpi=300, bbox_inches='tight')
        print("图片已保存为 'permutation_test_results.png'")
    
    plt.show()

def calculate_ref_base_ratios_with_permutation():
    """
    计算两个位点的ref base ratio并进行置换检验
    """
    
    # chr10:541528 数据 (A→G SNP, ref=A)
    chr10_data = {
        'position': 'chr10:541528',
        'ref_allele': 'A',
        'alt_allele': 'G',
        'AF_total': 0.309182,
        'AF_AFR': 0.657895,
        'AF_AMR': 0.162245,
        'AF_EAS': 0.267521,
        'AF_EUR': 0.1406,
        'AF_SAS': 0.128952
    }
    
    # chr1:1379664 数据 (G→A SNP, ref=G)
    chr1_data = {
        'position': 'chr1:1379664',
        'ref_allele': 'G',
        'alt_allele': 'A',
        'AF_total': 0.724547,
        'AF_AFR': 0.357783,
        'AF_AMR': 0.842857,
        'AF_EAS': 0.880342,
        'AF_EUR': 0.913112,
        'AF_SAS': 0.822795
    }
    
    # 计算ref base ratio (1 - AF)
    chr10_ref_ratios = {
        'total': 1 - chr10_data['AF_total'],
        'AFR': 1 - chr10_data['AF_AFR'],
        'AMR': 1 - chr10_data['AF_AMR'],
        'EAS': 1 - chr10_data['AF_EAS'],
        'EUR': 1 - chr10_data['AF_EUR'],
        'SAS': 1 - chr10_data['AF_SAS']
    }
    
    chr1_ref_ratios = {
        'total': 1 - chr1_data['AF_total'],
        'AFR': 1 - chr1_data['AF_AFR'],
        'AMR': 1 - chr1_data['AF_AMR'],
        'EAS': 1 - chr1_data['AF_EAS'],
        'EUR': 1 - chr1_data['AF_EUR'],
        'SAS': 1 - chr1_data['AF_SAS']
    }
    
    # 打印基本信息
    print("=" * 70)
    print("SNP位点信息:")
    print(f"位点1: {chr10_data['position']} ({chr10_data['ref_allele']}→{chr10_data['alt_allele']})")
    print(f"位点2: {chr1_data['position']} ({chr1_data['ref_allele']}→{chr1_data['alt_allele']})")
    print("=" * 70)
    
    # 打印ref base ratio
    print("\nRef Base Ratio (1 - AF):")
    print("-" * 50)
    populations = ['total', 'AFR', 'AMR', 'EAS', 'EUR', 'SAS']
    
    print(f"{'Population':<10} {'chr10:541528':<15} {'chr1:1379664':<15}")
    print("-" * 50)
    for pop in populations:
        print(f"{pop:<10} {chr10_ref_ratios[pop]:<15.6f} {chr1_ref_ratios[pop]:<15.6f}")
    
    # 准备计算相关系数的数据（使用5个人群）
    populations_for_corr = ['AFR', 'AMR', 'EAS', 'EUR', 'SAS']
    
    chr10_values = np.array([chr10_ref_ratios[pop] for pop in populations_for_corr])
    chr1_values = np.array([chr1_ref_ratios[pop] for pop in populations_for_corr])
    
    # 标准皮尔森相关系数
    correlation, p_value_pearson = pearsonr(chr10_values, chr1_values)
    
    print("\n" + "=" * 70)
    print("标准皮尔森相关系数:")
    print("=" * 70)
    print(f"相关系数 (r): {correlation:.6f}")
    print(f"标准P值: {p_value_pearson:.6f}")
    
    # 置换检验
    print("\n" + "=" * 70)
    print("置换检验:")
    print("=" * 70)
    
    observed_corr, perm_corrs, p_value_perm = permutation_test(
        chr10_values, chr1_values, n_permutations=10000)
    
    print(f"\n置换检验结果:")
    print(f"观察相关系数: {observed_corr:.6f}")
    print(f"置换检验P值: {p_value_perm:.6f}")
    
    # 计算效应量
    effect_size = abs(observed_corr)
    if effect_size >= 0.8:
        effect_desc = "大效应"
    elif effect_size >= 0.5:
        effect_desc = "中等效应"
    elif effect_size >= 0.3:
        effect_desc = "小效应"
    else:
        effect_desc = "微小效应"
    
    print(f"效应量: {effect_desc} (|r| = {effect_size:.3f})")
    
    # 置换检验的统计摘要
    print(f"\n置换分布统计:")
    print(f"均值: {np.mean(perm_corrs):.6f}")
    print(f"标准差: {np.std(perm_corrs):.6f}")
    print(f"95%置信区间: [{np.percentile(perm_corrs, 2.5):.4f}, {np.percentile(perm_corrs, 97.5):.4f}]")
    
    # 比较两种检验方法
    print("\n" + "=" * 70)
    print("检验方法比较:")
    print("=" * 70)
    print(f"{'方法':<15} {'P值':<12} {'显著性':<10}")
    print("-" * 40)
    print(f"{'Pearson检验':<15} {p_value_pearson:<12.6f} {'是' if p_value_pearson < 0.05 else '否':<10}")
    print(f"{'置换检验':<15} {p_value_perm:<12.6f} {'是' if p_value_perm < 0.05 else '否':<10}")
    
    # 最终结论
    print("\n" + "=" * 70)
    print("最终结论:")
    print("=" * 70)
    
    if p_value_perm < 0.001:
        significance = "极显著 (p < 0.001)"
    elif p_value_perm < 0.01:
        significance = "高度显著 (p < 0.01)"
    elif p_value_perm < 0.05:
        significance = "显著 (p < 0.05)"
    else:
        significance = "不显著 (p ≥ 0.05)"
    
    direction = "正相关" if observed_corr > 0 else "负相关"
    
    print(f"两个SNP位点的ref base ratio在不同人群中呈现{direction}关系")
    print(f"相关系数: r = {observed_corr:.4f}")
    print(f"统计显著性: {significance}")
    print(f"效应量: {effect_desc}")
    
    if abs(p_value_pearson - p_value_perm) > 0.01:
        print(f"\n注意: 标准检验和置换检验的P值存在差异，建议以置换检验结果为准")
    
    # 绘制结果图
    try:
        plot_permutation_results(observed_corr, perm_corrs, p_value_perm, save_plot=True)
    except ImportError:
        print("\n注意: 未安装matplotlib，无法绘制图形")
    except Exception as e:
        print(f"\n绘图时出现错误: {e}")
    
    return observed_corr, p_value_perm, perm_corrs

def main():
    """
    主函数
    """
    try:
        correlation, p_value_perm, perm_corrs = calculate_ref_base_ratios_with_permutation()
        return 0
        
    except Exception as e:
        print(f"分析过程中出现错误: {e}")
        return 1

if __name__ == "__main__":
    exit(main())