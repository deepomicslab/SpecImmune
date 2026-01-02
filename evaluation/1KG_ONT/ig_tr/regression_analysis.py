"""
Multivariate regression analysis for IG/TCR heterozygosity.

Analyzes heterozygous variant counts controlling for:
- Sequencing depth
- Read length
- Mapping quality
- Basecaller mode
- Population structure
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from matplotlib.backends.backend_pdf import PdfPages
import os

from vdj_1kgp_analysis import get_super_pop, get_sample_pop

OUTPUT_DIR = 'regression'
os.makedirs(OUTPUT_DIR, exist_ok=True)


def multivariate_regression_analysis(df, pdf):
    """
    Full regression model with all technical and population covariates.
    """
    print("\n" + "="*80)
    print("MULTIVARIATE REGRESSION ANALYSIS")
    print("="*80)
    
    # Prepare data - remove rows with missing values
    analysis_df = df.copy()
    
    # Create additional features
    analysis_df['log_depth'] = np.log10(analysis_df['depth'] + 1)
    
    # Use avg_read_length from the data if available, otherwise compute from allele lengths
    if 'avg_read_length' not in analysis_df.columns:
        analysis_df['avg_read_length'] = (analysis_df['length_1'] + analysis_df['length_2']) / 2
    
    # Create dummy variables for populations
    pop_dummies = pd.get_dummies(analysis_df['super_pop'], prefix='pop')
    analysis_df = pd.concat([analysis_df, pop_dummies], axis=1)
    
    # Create dummy variables for base_caller_mode if available
    if 'base_caller_mode' in analysis_df.columns:
        basecaller_dummies = pd.get_dummies(analysis_df['base_caller_mode'], prefix='basecaller')
        analysis_df = pd.concat([analysis_df, basecaller_dummies], axis=1)
        has_basecaller = True
    else:
        has_basecaller = False
    
    # Define predictors - include mapping_quality if available
    base_predictors = ['log_depth', 'avg_read_length']
    if 'mapping_quality' in analysis_df.columns:
        base_predictors.append('mapping_quality')
    
    # Remove rows with NaN values in the predictors
    required_cols = base_predictors + ['hete_variant_num']
    analysis_df = analysis_df.dropna(subset=required_cols)
    
    print(f"\nAnalyzing {len(analysis_df)} samples after removing missing values")
    print(f"Covariates included: {', '.join(base_predictors)}")
    if has_basecaller:
        print(f"Base caller modes: {analysis_df['base_caller_mode'].unique()}")
    
    pop_columns = [col for col in analysis_df.columns if col.startswith('pop_')]
    basecaller_columns = [col for col in analysis_df.columns if col.startswith('basecaller_')] if has_basecaller else []
    
    y = analysis_df['hete_variant_num']
    
    # Full Model: All covariates (technical + basecaller + population)
    print("\n" + "-"*80)
    print("Full Model: All Covariates (Technical + Base Caller + Population)")
    print("-"*80)
    all_predictors = base_predictors + pop_columns[:-1]  # Exclude one population to avoid multicollinearity
    if has_basecaller:
        all_predictors += basecaller_columns[:-1]
    X_full = analysis_df[all_predictors]
    model_full = LinearRegression()
    model_full.fit(X_full, y)
    r2_full = model_full.score(X_full, y)
    print(f"R² = {r2_full:.4f}")
    print(f"Number of predictors: {len(all_predictors)}")
    print("\nCoefficients:")
    for i, pred in enumerate(all_predictors):
        print(f"  {pred}: {model_full.coef_[i]:.4f}")
    
    # Model summary
    print("\n" + "="*80)
    print("MODEL SUMMARY")
    print("="*80)
    model_summary = pd.DataFrame({
        'Model': ['Full Model'],
        'R²': [r2_full],
        'Predictors': [f'{len(all_predictors)}']
    })
    print(model_summary)
    print(f"\nVariance explained: {r2_full*100:.2f}%")
    
    model_file = os.path.join(OUTPUT_DIR, 'model_summary.csv')
    model_summary.to_csv(model_file, index=False)
    
    # Analyze AFR vs non-AFR after controlling for technical covariates
    print("\n" + "="*80)
    print("AFR vs NON-AFR HETEROZYGOSITY COMPARISON")
    print("(After Controlling for Technical Factors)")
    print("="*80)
    
    # Create AFR indicator
    analysis_df['is_AFR'] = (analysis_df['super_pop'] == 'AFR').astype(int)
    
    # Predict using only technical covariates (no population)
    X_tech = analysis_df[base_predictors]
    if has_basecaller:
        X_tech = analysis_df[base_predictors + basecaller_columns[:-1]]
    model_tech = LinearRegression()
    model_tech.fit(X_tech, y)
    
    # Calculate residuals (heterozygosity after removing technical effects)
    analysis_df['residuals'] = y - model_tech.predict(X_tech)
    
    # Compare AFR vs non-AFR on residuals
    afr_residuals = analysis_df[analysis_df['is_AFR'] == 1]['residuals']
    non_afr_residuals = analysis_df[analysis_df['is_AFR'] == 0]['residuals']
    
    print(f"\nSample sizes:")
    print(f"  AFR: {len(afr_residuals)}")
    print(f"  Non-AFR: {len(non_afr_residuals)}")
    
    print(f"\nResiduals (heterozygosity after removing technical effects):")
    print(f"  AFR:     mean = {afr_residuals.mean():.4f} ± {afr_residuals.std():.4f}")
    print(f"  Non-AFR: mean = {non_afr_residuals.mean():.4f} ± {non_afr_residuals.std():.4f}")
    print(f"  Difference (AFR - Non-AFR): {afr_residuals.mean() - non_afr_residuals.mean():.4f}")
    
    # Statistical tests
    from scipy import stats
    t_stat, t_pval = stats.ttest_ind(afr_residuals, non_afr_residuals)
    u_stat, u_pval = stats.mannwhitneyu(afr_residuals, non_afr_residuals, alternative='two-sided')
    
    print(f"\nStatistical tests:")
    print(f"  T-test: t = {t_stat:.4f}, p = {t_pval:.2e}")
    print(f"  Mann-Whitney U: U = {u_stat:.0f}, p = {u_pval:.2e}")
    
    # Effect size (Cohen's d)
    pooled_std = np.sqrt(((len(afr_residuals)-1)*afr_residuals.std()**2 + 
                          (len(non_afr_residuals)-1)*non_afr_residuals.std()**2) / 
                         (len(afr_residuals) + len(non_afr_residuals) - 2))
    cohens_d = (afr_residuals.mean() - non_afr_residuals.mean()) / pooled_std
    print(f"  Cohen's d: {cohens_d:.4f}")
    
    # Conclusion
    print(f"\n{'='*80}")
    if t_pval < 0.001:
        sig_level = "p < 0.001"
    elif t_pval < 0.01:
        sig_level = "p < 0.01"
    elif t_pval < 0.05:
        sig_level = "p < 0.05"
    else:
        sig_level = f"p = {t_pval:.2f}"
    
    if afr_residuals.mean() > non_afr_residuals.mean() and t_pval < 0.05:
        conclusion = f"✓ AFR has SIGNIFICANTLY HIGHER heterozygosity than non-AFR ({sig_level})"
    elif afr_residuals.mean() < non_afr_residuals.mean() and t_pval < 0.05:
        conclusion = f"✗ AFR has SIGNIFICANTLY LOWER heterozygosity than non-AFR ({sig_level})"
    else:
        conclusion = f"○ No significant difference between AFR and non-AFR ({sig_level})"
    
    print(f"CONCLUSION: {conclusion}")
    print(f"{'='*80}")
    
    # Save comparison results
    comparison_df = pd.DataFrame({
        'group': ['AFR', 'Non-AFR'],
        'n_samples': [len(afr_residuals), len(non_afr_residuals)],
        'mean_residual': [afr_residuals.mean(), non_afr_residuals.mean()],
        'std_residual': [afr_residuals.std(), non_afr_residuals.std()],
        'median_residual': [afr_residuals.median(), non_afr_residuals.median()]
    })
    comparison_file = os.path.join(OUTPUT_DIR, 'afr_vs_nonafr_comparison.csv')
    comparison_df.to_csv(comparison_file, index=False)
    
    # Also show per-population breakdown
    print(f"\nPer-population residuals (for reference):")
    print("-" * 80)
    pop_residual_stats = []
    for pop in sorted(analysis_df['super_pop'].unique()):
        pop_residuals = analysis_df[analysis_df['super_pop'] == pop]['residuals']
        pop_residual_stats.append({
            'population': pop,
            'mean_residual': pop_residuals.mean(),
            'std_residual': pop_residuals.std(),
            'n_samples': len(pop_residuals)
        })
        print(f"  {pop}: mean = {pop_residuals.mean():.4f} ± {pop_residuals.std():.4f} (n={len(pop_residuals)})")
    
    residual_df = pd.DataFrame(pop_residual_stats)
    residual_df = residual_df.sort_values('mean_residual', ascending=False)
    residual_file = os.path.join(OUTPUT_DIR, 'population_residuals.csv')
    residual_df.to_csv(residual_file, index=False)
    
    # Plot: Two subplots comparing AFR vs non-AFR
    fig = plt.figure(figsize=(14, 5.5))
    gs = fig.add_gridspec(1, 2, hspace=0.25, wspace=0.3)
    
    # Plot 1: Violin plot (top left)
    ax1 = fig.add_subplot(gs[0, 0])
    
    # Remove extreme outliers for visualization (keep within 5th-95th percentile)
    non_afr_q5, non_afr_q95 = np.percentile(non_afr_residuals, [5, 95])
    afr_q5, afr_q95 = np.percentile(afr_residuals, [5, 95])
    
    non_afr_filtered = non_afr_residuals[(non_afr_residuals >= non_afr_q5) & (non_afr_residuals <= non_afr_q95)]
    afr_filtered = afr_residuals[(afr_residuals >= afr_q5) & (afr_residuals <= afr_q95)]
    
    parts = ax1.violinplot([non_afr_filtered, afr_filtered], positions=[1, 2], 
                           showmeans=False, showmedians=False, widths=0.65)
    
    # Style violin plots
    colors = ['#4A90E2', '#E85D5D']
    for i, pc in enumerate(parts['bodies']):
        pc.set_facecolor(colors[i])
        pc.set_alpha(0.8)
        pc.set_edgecolor('black')
        pc.set_linewidth(1.5)
    
    # Add mean markers and lines
    means = [non_afr_residuals.mean(), afr_residuals.mean()]
    ax1.scatter([1, 2], means, color='white', s=200, zorder=10, edgecolor='black', linewidth=2.5, marker='D')
    ax1.hlines(means, [0.7, 1.7], [1.3, 2.3], colors='black', linewidth=2.5, linestyle='-', zorder=9)
    
    # Labels and styling
    ax1.set_xticks([1, 2])
    ax1.set_xticklabels(['Non-AFR\n(n=283,170)', 'AFR\n(n=104,747)'], fontsize=12, fontweight='bold')
    ax1.axhline(y=0, color='gray', linestyle='--', linewidth=2, alpha=0.6)
    ax1.set_ylabel('Residual Heterozygosity', fontsize=13, fontweight='bold')
    ax1.set_title(f'AFR vs Non-AFR Comparison\np = {t_pval:.2e}', fontsize=14, fontweight='bold', pad=15)
    ax1.grid(True, alpha=0.3, axis='y', linewidth=0.5)
    ax1.set_xlim(0.5, 2.5)
    
    # Add mean value annotations
    for i, (pos, mean) in enumerate(zip([1, 2], means)):
        ax1.text(pos, mean + 0.2, f'{mean:.3f}', ha='center', va='bottom',
                fontsize=11, fontweight='bold', 
                bbox=dict(boxstyle='round,pad=0.4', facecolor='white', edgecolor='black', linewidth=1.5))
    
    # Plot 2: Histogram (top right)
    ax2 = fig.add_subplot(gs[0, 1])
    
    bins = np.linspace(min(non_afr_filtered.min(), afr_filtered.min()), 
                       max(non_afr_filtered.max(), afr_filtered.max()), 60)
    
    # Use step histograms for clarity
    ax2.hist(non_afr_filtered, bins=bins, alpha=0, histtype='step', 
             color=colors[0], linewidth=2.5, density=True)
    ax2.hist(afr_filtered, bins=bins, alpha=0, histtype='step',
             color=colors[1], linewidth=2.5, density=True)
    
    # Fill under curves with transparency
    counts_non, bins_non = np.histogram(non_afr_filtered, bins=bins, density=True)
    counts_afr, bins_afr = np.histogram(afr_filtered, bins=bins, density=True)
    ax2.fill_between(bins_non[:-1], counts_non, alpha=0.3, color=colors[0], step='post')
    ax2.fill_between(bins_afr[:-1], counts_afr, alpha=0.3, color=colors[1], step='post')
    
    # Add mean lines
    ax2.axvline(non_afr_residuals.mean(), color=colors[0], linestyle='--', linewidth=3, 
                label=f'Non-AFR: {non_afr_residuals.mean():.3f}', zorder=10)
    ax2.axvline(afr_residuals.mean(), color=colors[1], linestyle='--', linewidth=3, 
                label=f'AFR: {afr_residuals.mean():.3f}', zorder=10)
    ax2.axvline(0, color='gray', linestyle='-', linewidth=2, alpha=0.5)
    
    ax2.set_xlabel('Residual Heterozygosity', fontsize=13, fontweight='bold')
    ax2.set_ylabel('Density', fontsize=13, fontweight='bold')
    ax2.set_title(f'Distribution Comparison\nΔ = {afr_residuals.mean() - non_afr_residuals.mean():.4f}', 
                  fontsize=14, fontweight='bold', pad=15)
    ax2.legend(fontsize=10, loc='upper right', frameon=True, shadow=True)
    ax2.grid(True, alpha=0.3, axis='y', linewidth=0.5)
    
    plt.tight_layout()
    pdf.savefig(fig, dpi=150, bbox_inches='tight')
    plt.close()
    
    return analysis_df, model_full, model_summary

if __name__ == "__main__":
    super_pop_file = "../hla/20131219.populations.tsv"
    sample_pop_file = "../hla/20130606_sample_info.xlsx"
    allele_file = "1kg_merged_samples.ig_tr.new.csv"

    print("="*80)
    print("IG/TCR HETEROZYGOSITY REGRESSION ANALYSIS")
    print("Controlling for Depth and Technical Covariates")
    print("="*80)

    # Load data
    print("\nLoading data...")
    super_pop_dict = get_super_pop(super_pop_file)
    sample_pop_dict = get_sample_pop(sample_pop_file)

    df = pd.read_csv(allele_file)
    
    # Add population annotations
    pop_list = []
    super_pop_list = []
    for index, row in df.iterrows():
        if row['Sample'] not in sample_pop_dict:
            pop = 'Non-pop'
            super_pop = 'Non-pop'
        else:
            pop = sample_pop_dict[row['Sample']]
            super_pop = super_pop_dict[pop]
        pop_list.append(pop)
        super_pop_list.append(super_pop)
    
    df['pop'] = pop_list
    df['super_pop'] = super_pop_list
    
    # Remove non-population samples
    df = df[df['super_pop'] != 'Non-pop'].copy()
    
    print(f"\nDataset summary:")
    print(f"  Total samples: {len(df)}")
    print(f"  Populations: {df['super_pop'].unique()}")
    print(f"  Genes: {len(df['gene'].unique())}")
    
    # Create single PDF file for all plots
    output_pdf = os.path.join(OUTPUT_DIR, 'regression_analysis_complete.pdf')
    
    print("\n" + "="*80)
    print("RUNNING MULTIVARIATE REGRESSION ANALYSIS")
    print("="*80)
    print(f"\nAll plots will be saved to: {output_pdf}")
    
    with PdfPages(output_pdf) as pdf:
        # Multivariate regression analysis
        analysis_df, model_full, model_summary = multivariate_regression_analysis(df, pdf)
    
    print("\n" + "="*80)
    print("ANALYSIS COMPLETE")
    print("="*80)
    print(f"\n✓ All plots and analysis saved to: {output_pdf}")
    print(f"✓ Additional CSV files saved to: {OUTPUT_DIR}/")
    print("\n" + "="*80)
    print("RESULTS SUMMARY")
    print("="*80)
    print(f"\nFull Model Performance:")
    for idx, row in model_summary.iterrows():
        print(f"  R² = {row['R²']:.4f} with {row['Predictors']} predictors")
        print(f"  Variance explained: {row['R²']*100:.2f}%")
    print("\n" + "="*80)

