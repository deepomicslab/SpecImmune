import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from itertools import combinations
from matplotlib.ticker import FuncFormatter, MaxNLocator

# 设置全局样式
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 20
plt.rcParams['axes.linewidth'] = 1.2
plt.rcParams['xtick.major.width'] = 1.2
plt.rcParams['ytick.major.width'] = 1.2

# 读取数据
df_typing = pd.read_csv('typing_success_rates.csv')
df_acc_common = pd.read_csv('accuracy_common_samples.csv')
df_acc_penalty = pd.read_csv('accuracy_with_penalty.csv')
df_stats = pd.read_csv('statistical_comparisons.csv')

# ==================== 统一软件名称 ====================
software_name_mapping = {
    'SpecLong': 'Specimmune',
    'speclong': 'Specimmune',
    'SPECLONG': 'Specimmune',
    'SpecHLA': 'Spechla',
    'spechla': 'Spechla',
    'SPECHLA': 'Spechla',
    'HLA*LA': 'HLA*LA',
    'HLALA': 'HLA*LA',
    'HLA-LA': 'HLA*LA'
}

# 对有 Software 列的数据框应用名称映射
for df in [df_typing, df_acc_common, df_acc_penalty]:
    if 'Software' in df.columns:
        df['Software'] = df['Software'].str.strip()
        df['Software'] = df['Software'].replace(software_name_mapping)

# 对 df_stats 的 Comparison 列进行处理
if 'Comparison' in df_stats.columns:
    df_stats['Comparison'] = df_stats['Comparison'].str.strip()
    for old_name, new_name in software_name_mapping.items():
        df_stats['Comparison'] = df_stats['Comparison'].str.replace(old_name, new_name, regex=False)

# 软件颜色映射
SOFTWARE_COLORS = {
    "HLA*LA": "#98B6C4",
    "Spechla": "#8D93AF", 
    "Specimmune": "#3E4271"
}

SOFTWARE_ORDER = ["HLA*LA", "Spechla", "Specimmune"]

# ==================== 数据集名称映射和顺序 ====================
# 数据中的实际名称 -> 显示名称
DATASET_DISPLAY_NAMES = {
    'hprc_hifi': 'HPRC HiFi',
    'hprc_ont': 'HPRC ONT',
    'hgsvc2_hifi': 'HGSVC HiFi',
    'hgsvc2_clr': 'HGSVC CLR',
    '1kg_ont': '1kGP ONT'  # 新增
}

# 按照您要求的顺序 (新增 1kg_ont)
DATASET_ORDER = ['hprc_hifi', 'hprc_ont', 'hgsvc2_hifi', 'hgsvc2_clr', '1kg_ont']

print("\n" + "="*80)
print("数据预处理完成 - 软件名称已统一")
print("="*80)
print(f"统一后的软件名称: {SOFTWARE_ORDER}")
print(f"数据集显示顺序: {[DATASET_DISPLAY_NAMES[ds] for ds in DATASET_ORDER]}")
print(f"df_typing 中的软件: {sorted(df_typing['Software'].unique())}")
print(f"df_acc_common 中的软件: {sorted(df_acc_common['Software'].unique())}")
if 'Dataset' in df_acc_common.columns:
    print(f"df_acc_common 中的数据集: {sorted(df_acc_common['Dataset'].unique())}")
print("="*80 + "\n")

# ==================== 辅助函数：添加显著性标注 ====================
def add_significance_brackets(ax, x_positions, values, comparisons_data, y_offset=0.02, bracket_height=0.03):
    """在图表上添加显著性标注"""
    
    valid_values = [v for v in values if not np.isnan(v)]
    if len(valid_values) == 0:
        return
    
    y_max = max(valid_values)
    
    pairs = list(combinations(range(len(SOFTWARE_ORDER)), 2))
    pairs_sorted = sorted(pairs, key=lambda p: abs(p[1] - p[0]))
    
    level = 0
    for i, j in pairs_sorted:
        software1 = SOFTWARE_ORDER[i]
        software2 = SOFTWARE_ORDER[j]
        
        comparison_name = f"{software1} vs {software2}"
        comparison_reverse = f"{software2} vs {software1}"
        
        result = comparisons_data[
            (comparisons_data['Comparison'] == comparison_name) |
            (comparisons_data['Comparison'] == comparison_reverse)
        ]
        
        if len(result) == 0:
            continue
        
        is_significant = result['Fisher_Significant'].values[0]
        p_value = result['Fisher_P'].values[0]
        
        if is_significant:
            x1, x2 = x_positions[i], x_positions[j]
            y = y_max + y_offset + level * bracket_height
            
            ax.plot([x1, x2], [y, y], 'k-', linewidth=1.2)
            ax.plot([x1, x1], [y - 0.01, y], 'k-', linewidth=1.2)
            ax.plot([x2, x2], [y - 0.01, y], 'k-', linewidth=1.2)
            
            if p_value < 0.001:
                sig_text = '***'
            elif p_value < 0.01:
                sig_text = '**'
            elif p_value < 0.05:
                sig_text = '*'
            else:
                sig_text = 'ns'
            
            ax.text((x1 + x2) / 2, y + 0.005, sig_text, 
                   ha='center', va='bottom', fontsize=10, fontweight='bold')
            
            level += 1

# ==================== 图1: 分型成功率箱线图 ====================
def plot_typing_success_rate():
    """绘制分型成功率箱线图"""
    fig, axes = plt.subplots(1, 3, figsize=(18, 6), dpi=300)
    
    class_names = ['Class I Genes', 'Class II Genes', 'All Genes']
    
    for idx, class_name in enumerate(class_names):
        ax = axes[idx]
        
        plot_data = df_typing[df_typing['Class'] == class_name].copy()
        
        if len(plot_data) == 0:
            continue
        
        plot_data['Software'] = pd.Categorical(
            plot_data['Software'], 
            categories=SOFTWARE_ORDER, 
            ordered=True
        )
        plot_data = plot_data.sort_values('Software')
        
        sns.boxplot(
            data=plot_data,
            x='Software',
            y='Typing_Success_Rate',
            ax=ax,
            order=SOFTWARE_ORDER,
            hue='Software',
            palette=SOFTWARE_COLORS,
            legend=False,
            showmeans=True,
            meanprops=dict(marker='D', markerfacecolor='red', markersize=6, zorder=3)
        )
        
        sns.stripplot(
            data=plot_data,
            x='Software',
            y='Typing_Success_Rate',
            ax=ax,
            order=SOFTWARE_ORDER,
            color='black',
            alpha=0.4,
            size=5,
            jitter=0.2
        )
        
        means = []
        for software in SOFTWARE_ORDER:
            software_data = plot_data[plot_data['Software'] == software]
            if len(software_data) > 0:
                means.append(software_data['Typing_Success_Rate'].mean())
            else:
                means.append(np.nan)
        
        stat_data = df_stats[df_stats['Class'] == class_name]
        if len(stat_data) > 0:
            add_significance_brackets(
                ax, 
                x_positions=list(range(len(SOFTWARE_ORDER))),
                values=means,
                comparisons_data=stat_data,
                y_offset=0.05,
                bracket_height=0.05
            )
        
        ax.set_title(f'{class_name}\nTyping Success Rate', fontsize=12, fontweight='bold')
        ax.set_xlabel('Software', fontsize=15)
        ax.set_ylabel('Typing Success Rate', fontsize=11)
        ax.set_ylim([0, 1.25])
        
        ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: '{:.0%}'.format(y)))
        ax.tick_params(axis='x', rotation=45)
        
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.yaxis.grid(True, alpha=0.3, linestyle='--')
        ax.set_axisbelow(True)
    
    plt.tight_layout()
    plt.savefig('typing_success_rate_boxplot.png', dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig('typing_success_rate_boxplot.pdf', bbox_inches='tight', facecolor='white')
    print("✓ 已保存: typing_success_rate_boxplot.png/pdf")
    plt.close()

# ==================== 图2: 准确率箱线图（共同样本） ====================
def plot_accuracy_common():
    """绘制共同样本准确率箱线图"""
    fig, axes = plt.subplots(1, 3, figsize=(18, 6), dpi=300)
    
    class_names = ['Class I Genes', 'Class II Genes', 'All Genes']
    
    for idx, class_name in enumerate(class_names):
        ax = axes[idx]
        
        plot_data = df_acc_common[df_acc_common['Class'] == class_name].copy()
        
        if len(plot_data) == 0:
            continue
        
        plot_data['Software'] = pd.Categorical(
            plot_data['Software'], 
            categories=SOFTWARE_ORDER, 
            ordered=True
        )
        plot_data = plot_data.sort_values('Software')
        
        sns.boxplot(
            data=plot_data,
            x='Software',
            y='Accuracy_Common',
            ax=ax,
            order=SOFTWARE_ORDER,
            hue='Software',
            palette=SOFTWARE_COLORS,
            legend=False,
            showmeans=True,
            meanprops=dict(marker='D', markerfacecolor='red', markersize=6, zorder=3)
        )
        
        sns.stripplot(
            data=plot_data,
            x='Software',
            y='Accuracy_Common',
            ax=ax,
            order=SOFTWARE_ORDER,
            color='black',
            alpha=0.4,
            size=5,
            jitter=0.2
        )
        
        means = []
        for software in SOFTWARE_ORDER:
            software_data = plot_data[plot_data['Software'] == software]
            if len(software_data) > 0:
                means.append(software_data['Accuracy_Common'].mean())
            else:
                means.append(np.nan)
        
        stat_data = df_stats[df_stats['Class'] == class_name]
        if len(stat_data) > 0:
            add_significance_brackets(
                ax, 
                x_positions=list(range(len(SOFTWARE_ORDER))),
                values=means,
                comparisons_data=stat_data,
                y_offset=0.03,
                bracket_height=0.04
            )
        
        ax.set_title(f'{class_name}\nAccuracy (Common Samples)', fontsize=12, fontweight='bold')
        ax.set_xlabel('Software', fontsize=11)
        ax.set_ylabel('Accuracy', fontsize=11)
        ax.set_ylim([0.5, 1.2])
        
        ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: '{:.0%}'.format(y)))
        ax.tick_params(axis='x', rotation=45)
        
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.yaxis.grid(True, alpha=0.3, linestyle='--')
        ax.set_axisbelow(True)
    
    plt.tight_layout()
    plt.savefig('accuracy_common_boxplot.png', dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig('accuracy_common_boxplot.pdf', bbox_inches='tight', facecolor='white')
    print("✓ 已保存: accuracy_common_boxplot.png/pdf")
    plt.close()

# ==================== 图3: 准确率箱线图（含惩罚） ====================
def plot_accuracy_penalty():
    """绘制含惩罚准确率箱线图"""
    fig, axes = plt.subplots(1, 3, figsize=(18, 6), dpi=300)
    
    class_names = ['Class I Genes', 'Class II Genes', 'All Genes']
    
    for idx, class_name in enumerate(class_names):
        ax = axes[idx]
        
        plot_data = df_acc_penalty[df_acc_penalty['Class'] == class_name].copy()
        
        if len(plot_data) == 0:
            continue
        
        plot_data['Software'] = pd.Categorical(
            plot_data['Software'], 
            categories=SOFTWARE_ORDER, 
            ordered=True
        )
        plot_data = plot_data.sort_values('Software')
        
        sns.boxplot(
            data=plot_data,
            x='Software',
            y='Accuracy_With_Penalty',
            ax=ax,
            order=SOFTWARE_ORDER,
            hue='Software',
            palette=SOFTWARE_COLORS,
            legend=False,
            showmeans=True,
            meanprops=dict(marker='D', markerfacecolor='red', markersize=6, zorder=3)
        )
        
        sns.stripplot(
            data=plot_data,
            x='Software',
            y='Accuracy_With_Penalty',
            ax=ax,
            order=SOFTWARE_ORDER,
            color='black',
            alpha=0.4,
            size=5,
            jitter=0.2
        )
        
        means = []
        for software in SOFTWARE_ORDER:
            software_data = plot_data[plot_data['Software'] == software]
            if len(software_data) > 0:
                means.append(software_data['Accuracy_With_Penalty'].mean())
            else:
                means.append(np.nan)
        
        stat_data = df_stats[df_stats['Class'] == class_name]
        if len(stat_data) > 0:
            add_significance_brackets(
                ax, 
                x_positions=list(range(len(SOFTWARE_ORDER))),
                values=means,
                comparisons_data=stat_data,
                y_offset=0.03,
                bracket_height=0.04
            )
        
        ax.set_title(f'{class_name}\nAccuracy (With Typing Penalty)', fontsize=12, fontweight='bold')
        ax.set_xlabel('Software', fontsize=11)
        ax.set_ylabel('Accuracy', fontsize=11)
        ax.set_ylim([0.5, 1.2])
        
        ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: '{:.0%}'.format(y)))
        ax.tick_params(axis='x', rotation=45)
        
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.yaxis.grid(True, alpha=0.3, linestyle='--')
        ax.set_axisbelow(True)
    
    plt.tight_layout()
    plt.savefig('accuracy_penalty_boxplot.png', dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig('accuracy_penalty_boxplot.pdf', bbox_inches='tight', facecolor='white')
    print("✓ 已保存: accuracy_penalty_boxplot.png/pdf")
    plt.close()

# ==================== 图4: 按数据集分组（扩展为2x3布局支持5个数据集） ====================
def plot_by_dataset():
    """按数据集分组展示准确率 - 2x3布局支持5个数据集"""
    
    datasets = DATASET_ORDER
    
    # 改为 2 行 3 列布局
    fig, axes = plt.subplots(2, 3, figsize=(12, 8), dpi=300)
    axes = axes.flatten()
    
    for idx, dataset in enumerate(datasets):
        if idx >= 6:  # 最多6个子图
            break
            
        ax = axes[idx]
        
        dataset_data = df_acc_common[
            (df_acc_common['Dataset'] == dataset) &
            (df_acc_common['Class'].isin(['Class I Genes', 'Class II Genes']))
        ].copy()
        
        if len(dataset_data) == 0:
            ax.text(0.5, 0.5, f'No data for\n{DATASET_DISPLAY_NAMES.get(dataset, dataset)}', 
                   ha='center', va='center', fontsize=15, transform=ax.transAxes)
            ax.set_xticks([])
            ax.set_yticks([])
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['left'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            continue
        
        class_labels = ['Class I', 'Class II']
        x = np.array([0, 1.02])
        width = 0.25
        
        all_bars_data = {}
        
        for i, software in enumerate(SOFTWARE_ORDER):
            software_data = dataset_data[dataset_data['Software'] == software]
            
            if len(software_data) == 0:
                continue
            
            accuracies = []
            errors_low = []
            errors_upp = []
            
            for class_name in ['Class I Genes', 'Class II Genes']:
                class_row = software_data[software_data['Class'] == class_name]
                if len(class_row) > 0:
                    acc = class_row['Accuracy_Common'].values[0]
                    ci_low = class_row['CI_Low'].values[0]
                    ci_upp = class_row['CI_Upp'].values[0]
                    accuracies.append(acc)
                    errors_low.append(acc - ci_low)
                    errors_upp.append(ci_upp - acc)
                else:
                    accuracies.append(np.nan)
                    errors_low.append(0)
                    errors_upp.append(0)
            
            valid_indices = [j for j, acc in enumerate(accuracies) if not np.isnan(acc)]
            
            if len(valid_indices) > 0:
                bars = ax.bar(
                    x[valid_indices] + (i - 1) * width, 
                    [accuracies[j] for j in valid_indices], 
                    width,
                    label=software if idx == 0 else "",
                    color=SOFTWARE_COLORS[software],
                    edgecolor='white',
                    linewidth=1.5,
                    yerr=[[errors_low[j] for j in valid_indices], 
                        [errors_upp[j] for j in valid_indices]],
                    capsize=4,
                    error_kw={
                        'elinewidth': 1.2,
                        'capthick': 1.2,
                        'ecolor': 'black',
                        'alpha': 0.6
                    }
                )
                
                # 添加样本数量标注
                for bar_idx, j in enumerate(valid_indices):
                    class_name_full = ['Class I Genes', 'Class II Genes'][j]
                    class_row = software_data[software_data['Class'] == class_name_full]
                    
                    if len(class_row) > 0 and 'Common_Samples' in class_row.columns:
                        n_value = class_row['Common_Samples'].values[0]
                        bar = bars[bar_idx]
                        height = bar.get_height()
                        error_upper = errors_upp[j]
                        
                        ax.text(
                            bar.get_x() + bar.get_width() / 2,
                            height + error_upper + 0.015,
                            f'n={int(n_value)}',
                            ha='center',
                            va='bottom',
                            fontsize=10,
                            # color='dimgray',
                            style='italic'
                        )
                
                all_bars_data[software] = {
                    'positions': x + (i - 1) * width,
                    'heights': accuracies,
                    'valid_indices': valid_indices
                }
                    
        # 添加显著性标注
        for class_idx, class_name in enumerate(['Class I Genes', 'Class II Genes']):
            class_stat_data = df_stats[
                (df_stats['Class'] == class_name) & 
                (df_stats['Dataset'] == dataset)
            ]
            
            if len(class_stat_data) == 0:
                continue
            
            heights_for_class = []
            positions_for_class = []
            software_indices = []
            
            for sw_idx, software in enumerate(SOFTWARE_ORDER):
                if software in all_bars_data:
                    if class_idx in all_bars_data[software]['valid_indices']:
                        positions_for_class.append(all_bars_data[software]['positions'][class_idx])
                        heights_for_class.append(all_bars_data[software]['heights'][class_idx])
                        software_indices.append(sw_idx)
            
            if len(heights_for_class) >= 2:
                y_max = max([h for h in heights_for_class if not np.isnan(h)])+0.02
                
                pairs = list(combinations(range(len(software_indices)), 2))
                level = 0
                
                for i, j in pairs:
                    actual_i = software_indices[i]
                    actual_j = software_indices[j]
                    
                    software1 = SOFTWARE_ORDER[actual_i]
                    software2 = SOFTWARE_ORDER[actual_j]
                    
                    comparison_name = f"{software1} vs {software2}"
                    comparison_reverse = f"{software2} vs {software1}"
                    
                    result = class_stat_data[
                        (class_stat_data['Comparison'] == comparison_name) |
                        (class_stat_data['Comparison'] == comparison_reverse)
                    ]
                    
                    if len(result) > 0 and result['Fisher_Significant'].values[0]:
                        p_value = result['Fisher_P'].values[0]
                        
                        x1, x2 = positions_for_class[i], positions_for_class[j]
                        y = y_max + 0.03 + level * 0.04
                        
                        ax.plot([x1, x2], [y, y], 'k-', linewidth=1.2)
                        ax.plot([x1, x1], [y - 0.01, y], 'k-', linewidth=1.2)
                        ax.plot([x2, x2], [y - 0.01, y], 'k-', linewidth=1.2)
                        
                        if p_value < 0.001:
                            sig_text = '***'
                        elif p_value < 0.01:
                            sig_text = '**'
                        elif p_value < 0.05:
                            sig_text = '*'
                        else:
                            sig_text = 'ns'
                        
                        ax.text((x1 + x2) / 2, y + 0.005, sig_text, 
                               ha='center', va='bottom', fontsize=10, fontweight='bold')
                        
                        level += 1
        
        # 使用显示名称作为标题
        display_name = DATASET_DISPLAY_NAMES.get(dataset, dataset)
        ax.set_title(display_name, fontsize=15, fontweight='bold', pad=10)
        ax.set_ylim([0.5, 1.15])
        ax.set_ylabel('Accuracy', fontsize=15, fontweight='bold')
        
        # ========== Y轴优化：只显示 ≤100% 的刻度 ==========
        ax.yaxis.set_major_locator(MaxNLocator(nbins=6, prune='upper'))
        
        def format_yticks(y, pos):
            if y <= 1.01:
                return f'{y:.0%}'
            else:
                return ''
        
        ax.yaxis.set_major_formatter(FuncFormatter(format_yticks))
        # ================================================
        
        ax.set_xticks(x)
        ax.set_xticklabels(class_labels, fontsize=12, fontweight='bold')
        ax.set_xlabel('')
        
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.yaxis.grid(True, linestyle='--', alpha=0.3, linewidth=1)
        ax.set_axisbelow(True)
        
        if idx == 0:
            ax.legend(
                frameon=False, 
                loc='lower right',
                fontsize=12,
                handlelength=1.5,
                handleheight=1.2
            )
    
    # 隐藏多余的子图（如果有）
    for idx in range(len(datasets), 6):
        axes[idx].axis('off')
    
    plt.tight_layout()
    plt.savefig('accuracy_by_dataset.png', dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig('accuracy_by_dataset.pdf', bbox_inches='tight', facecolor='white')
    print("✓ 已保存: accuracy_by_dataset.png/pdf")
    plt.close()

# ==================== 主函数 ====================
def main():
    """生成所有可视化"""
    print("\n" + "="*80)
    print("开始生成高质量可视化图表（带显著性标注）...")
    print("="*80 + "\n")
    
    plot_typing_success_rate()
    plot_accuracy_common()
    plot_accuracy_penalty()
    plot_by_dataset()
    
    print("\n" + "="*80)
    print("主要图表生成完成！")
    print("="*80)
    print("\n生成的文件：")
    print("1. typing_success_rate_boxplot.png/pdf - 分型成功率箱线图（带显著性）")
    print("2. accuracy_common_boxplot.png/pdf - 共同样本准确率箱线图（带显著性）")
    print("3. accuracy_penalty_boxplot.png/pdf - 含惩罚准确率箱线图（带显著性）")
    print("4. accuracy_by_dataset.png/pdf - 按数据集分组的条形图（2x3布局）")
    print(f"\n   数据集显示顺序: {' → '.join([DATASET_DISPLAY_NAMES[ds] for ds in DATASET_ORDER])}")
    print(f"   总共 {len(DATASET_ORDER)} 个数据集")
    print("\n配色方案：")
    print("  HLA*LA: #98B6C4 (浅蓝)")
    print("  Spechla: #8D93AF (紫灰)")
    print("  Specimmune: #3E4271 (深紫)")
    print("\n✓ Y轴优化：大于100%的刻度标签已隐藏")

if __name__ == "__main__":
    main()