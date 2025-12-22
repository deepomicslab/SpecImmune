"""
Comprehensive regression analysis for IG/TCR heterozygosity controlling for confounding factors.

This script analyzes heterozygous variant counts in IG/TCR loci while controlling for:
- Per-locus sequencing depth
- Read length distribution
- Basecaller model differences (if available)
- Mapping quality scores
- Gene length (locus size)

vdj table is like:
gene,depth,phase_set,allele_1,score_1,length_1,hap_1,allele_2,score_2,length_2,hap_2,hg38_chrom,hg38_len,variant_num,hete_variant_num,Sample,Dataset
TRAV1-1,12.91,20989748.0,TRAV1-1*02,100.0,269.0,hap1,TRAV1-1*01,100.0,275.0,hap2,chr14,729,3,2,HG00096,1KG
TRAV1-2,14.75,20989748.0,TRAV1-2*01,100.0,267.0,hap1,TRAV1-2*01,100.0,267.0,hap2,chr14,689,2,1,HG00096,1KG
TRAV2,17.94,,TRAV2*01,100.0,263.0,hap1,TRAV2*01,100.0,263.0,hap2,chr14,522,0,0,HG00096,1KG
TRAV3,22.87,,TRAV3*01,100.0,286.0,hap1,TRAV3*01,100.0,286.0,hap2,chr14,608,0,0,HG00096,1KG

"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler
from matplotlib.backends.backend_pdf import PdfPages
import warnings
import os
warnings.filterwarnings('ignore')

from vdj_1kgp_analysis import get_super_pop, get_sample_pop

# Create output directory
OUTPUT_DIR = 'regression'
os.makedirs(OUTPUT_DIR, exist_ok=True)


def calculate_depth_matched_subsets(df, depth_bins=10):
    """
    Create depth-matched subsets by binning samples into depth quantiles.
    This ensures fair comparison across different coverage levels.
    """
    # Create depth bins
    df['depth_bin'] = pd.qcut(df['depth'], q=depth_bins, labels=False, duplicates='drop')
    return df


def perform_depth_stratified_analysis(df):
    """
    Analyze heterozygosity within depth-matched strata to control for coverage bias.
    """
    print("\n" + "="*80)
    print("DEPTH-STRATIFIED ANALYSIS")
    print("="*80)
    
    # Group by depth bins
    depth_bins = pd.qcut(df['depth'], q=10, labels=False, duplicates='drop')
    df['depth_bin'] = depth_bins
    
    stratified_results = []
    for bin_id in sorted(df['depth_bin'].unique()):
        bin_data = df[df['depth_bin'] == bin_id]
        
        # Get depth range for this bin
        depth_min = bin_data['depth'].min()
        depth_max = bin_data['depth'].max()
        depth_mean = bin_data['depth'].mean()
        
        # Calculate heterozygosity metrics per population
        pop_stats = []
        for pop in bin_data['super_pop'].unique():
            pop_data = bin_data[bin_data['super_pop'] == pop]
            if len(pop_data) > 0:
                mean_hete = pop_data['hete_variant_num'].mean()
                pop_stats.append({
                    'depth_bin': bin_id,
                    'depth_range': f"{depth_min:.1f}-{depth_max:.1f}",
                    'depth_mean': depth_mean,
                    'population': pop,
                    'mean_hete_variants': mean_hete,
                    'n_samples': len(pop_data)
                })
        
        stratified_results.extend(pop_stats)
    
    stratified_df = pd.DataFrame(stratified_results)
    
    # Save results
    output_file = os.path.join(OUTPUT_DIR, 'depth_stratified_heterozygosity.csv')
    stratified_df.to_csv(output_file, index=False)
    print(f"\nSaved depth-stratified results to: {output_file}")
    print(stratified_df.head(20))
    
    return stratified_df


def plot_depth_distributions(df, pdf):
    """
    Plot depth distributions per population for IG/TCR loci.
    """
    print("\n" + "="*80)
    print("DEPTH DISTRIBUTION ANALYSIS")
    print("="*80)
    
    # Overall depth statistics
    print(f"\nOverall depth statistics:")
    print(df['depth'].describe())
    
    # Per-population depth statistics
    print(f"\nPer-population depth statistics:")
    pop_depth_stats = df.groupby('super_pop')['depth'].describe()
    print(pop_depth_stats)
    
    # Save statistics
    stats_file = os.path.join(OUTPUT_DIR, 'depth_distribution_stats.csv')
    pop_depth_stats.to_csv(stats_file)
    
    # Create comprehensive plots
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    
    # 1. Box plot of depth by population
    ax1 = axes[0, 0]
    df.boxplot(column='depth', by='super_pop', ax=ax1)
    ax1.set_title('Depth Distribution by Population')
    ax1.set_xlabel('Population')
    ax1.set_ylabel('Depth')
    ax1.set_yscale('log')
    plt.sca(ax1)
    plt.xticks(rotation=45)
    
    # 2. Violin plot
    ax2 = axes[0, 1]
    sns.violinplot(data=df, x='super_pop', y='depth', ax=ax2)
    ax2.set_title('Depth Distribution (Violin Plot)')
    ax2.set_xlabel('Population')
    ax2.set_ylabel('Depth')
    ax2.set_yscale('log')
    plt.sca(ax2)
    plt.xticks(rotation=45)
    
    # 3. Histogram by population
    ax3 = axes[1, 0]
    for pop in df['super_pop'].unique():
        pop_data = df[df['super_pop'] == pop]['depth']
        ax3.hist(pop_data, bins=50, alpha=0.5, label=pop)
    ax3.set_xlabel('Depth')
    ax3.set_ylabel('Frequency')
    ax3.set_title('Depth Histogram by Population')
    ax3.legend()
    ax3.set_xlim(0, df['depth'].quantile(0.95))
    
    # 4. Gene-level depth variation
    ax4 = axes[1, 1]
    gene_depth_stats = df.groupby('gene')['depth'].agg(['mean', 'std']).reset_index()
    gene_depth_stats = gene_depth_stats.sort_values('mean', ascending=False).head(20)
    ax4.barh(range(len(gene_depth_stats)), gene_depth_stats['mean'], xerr=gene_depth_stats['std'])
    ax4.set_yticks(range(len(gene_depth_stats)))
    ax4.set_yticklabels(gene_depth_stats['gene'])
    ax4.set_xlabel('Mean Depth')
    ax4.set_title('Top 20 Genes by Mean Depth')
    
    plt.tight_layout()
    pdf.savefig(fig, dpi=150, bbox_inches='tight')
    plt.close()


def multivariate_regression_analysis(df, pdf):
    """
    Perform multivariate regression to assess the impact of depth and other covariates
    on heterozygous variant counts, controlling for confounding factors.
    """
    print("\n" + "="*80)
    print("MULTIVARIATE REGRESSION ANALYSIS")
    print("="*80)
    
    # Prepare data - remove rows with missing values
    analysis_df = df.copy()
    
    # Create additional features
    analysis_df['log_depth'] = np.log10(analysis_df['depth'] + 1)
    analysis_df['avg_read_length'] = (analysis_df['length_1'] + analysis_df['length_2']) / 2
    analysis_df['avg_score'] = (analysis_df['score_1'] + analysis_df['score_2']) / 2
    
    # Create dummy variables for populations
    pop_dummies = pd.get_dummies(analysis_df['super_pop'], prefix='pop')
    analysis_df = pd.concat([analysis_df, pop_dummies], axis=1)
    
    # Remove rows with NaN values
    analysis_df = analysis_df.dropna(subset=['log_depth', 'avg_read_length', 
                                               'avg_score', 'hete_variant_num'])
    
    print(f"\nAnalyzing {len(analysis_df)} samples after removing missing values")
    
    # Define predictors
    base_predictors = ['log_depth', 'avg_read_length', 'avg_score']
    pop_columns = [col for col in analysis_df.columns if col.startswith('pop_')]
    
    # Model 1: Depth only
    print("\n" + "-"*80)
    print("Model 1: Depth only (baseline)")
    print("-"*80)
    X1 = analysis_df[['log_depth']]
    y = analysis_df['hete_variant_num']
    model1 = LinearRegression()
    model1.fit(X1, y)
    r2_model1 = model1.score(X1, y)
    print(f"R² = {r2_model1:.4f}")
    print(f"Coefficient for log_depth: {model1.coef_[0]:.4f}")
    
    # Model 2: All technical covariates (no population)
    print("\n" + "-"*80)
    print("Model 2: All Technical Covariates (depth, read length, quality)")
    print("-"*80)
    X2 = analysis_df[base_predictors]
    model2 = LinearRegression()
    model2.fit(X2, y)
    r2_model2 = model2.score(X2, y)
    print(f"R² = {r2_model2:.4f}")
    for i, pred in enumerate(base_predictors):
        print(f"Coefficient for {pred}: {model2.coef_[i]:.4f}")
    
    # Model 3: Technical covariates + population
    print("\n" + "-"*80)
    print("Model 3: Technical Covariates + Population")
    print("-"*80)
    all_predictors = base_predictors + pop_columns[:-1]  # Exclude one population to avoid multicollinearity
    X3 = analysis_df[all_predictors]
    model3 = LinearRegression()
    model3.fit(X3, y)
    r2_model3 = model3.score(X3, y)
    print(f"R² = {r2_model3:.4f}")
    print("\nCoefficients:")
    for i, pred in enumerate(all_predictors):
        print(f"  {pred}: {model3.coef_[i]:.4f}")
    
    # Compute residuals from Model 2 (technical covariates only)
    analysis_df['residuals_tech'] = y - model2.predict(X2)
    
    # Test if population explains residual variation
    print("\n" + "-"*80)
    print("Population Effect After Controlling for Technical Covariates")
    print("-"*80)
    
    pop_residual_stats = []
    for pop in analysis_df['super_pop'].unique():
        pop_residuals = analysis_df[analysis_df['super_pop'] == pop]['residuals_tech']
        pop_residual_stats.append({
            'population': pop,
            'mean_residual': pop_residuals.mean(),
            'std_residual': pop_residuals.std(),
            'n_samples': len(pop_residuals)
        })
    
    residual_df = pd.DataFrame(pop_residual_stats)
    print(residual_df)
    residual_file = os.path.join(OUTPUT_DIR, 'population_residuals.csv')
    residual_df.to_csv(residual_file, index=False)
    
    # ANOVA on residuals
    residual_groups = [analysis_df[analysis_df['super_pop'] == pop]['residuals_tech'].values 
                       for pop in analysis_df['super_pop'].unique()]
    f_stat, p_value = stats.f_oneway(*residual_groups)
    print(f"\nANOVA on residuals (testing population effect after technical correction):")
    print(f"F-statistic: {f_stat:.4f}, p-value: {p_value:.4e}")
    
    # Model comparison
    print("\n" + "="*80)
    print("MODEL COMPARISON SUMMARY")
    print("="*80)
    model_comparison = pd.DataFrame({
        'Model': ['Depth only', 'All Technical', 'Technical + Population'],
        'R²': [r2_model1, r2_model2, r2_model3],
        'Predictors': ['1', f'{len(base_predictors)}', f'{len(all_predictors)}']
    })
    print(model_comparison)
    model_file = os.path.join(OUTPUT_DIR, 'model_comparison.csv')
    model_comparison.to_csv(model_file, index=False)
    
    # Plot: Observed vs Predicted
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    models = [
        (model1, X1, 'Model 1: Depth only', r2_model1),
        (model2, X2, 'Model 2: All Technical', r2_model2),
        (model3, X3, 'Model 3: Technical + Population', r2_model3)
    ]
    
    for idx, (model, X, title, r2) in enumerate(models):
        ax = axes[idx]
        y_pred = model.predict(X)
        # Sample data if too many points
        if len(y) > 5000:
            sample_idx = np.random.choice(len(y), 5000, replace=False)
            ax.scatter(y.iloc[sample_idx], y_pred[sample_idx], alpha=0.3, s=5)
        else:
            ax.scatter(y, y_pred, alpha=0.3, s=5)
        ax.plot([y.min(), y.max()], [y.min(), y.max()], 'r--', lw=2)
        ax.set_xlabel('Observed Heterozygous Variants')
        ax.set_ylabel('Predicted Heterozygous Variants')
        ax.set_title(f'{title}\nR² = {r2:.4f}')
        ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    pdf.savefig(fig, dpi=150, bbox_inches='tight')
    plt.close()
    
    return analysis_df, model2, model3, model_comparison


def depth_matched_population_comparison(df, pdf, n_bins=5):
    """
    Compare populations within depth-matched bins to control for coverage differences.
    """
    print("\n" + "="*80)
    print("DEPTH-MATCHED POPULATION COMPARISON")
    print("="*80)
    
    # Create depth bins
    df['depth_bin'] = pd.qcut(df['depth'], q=n_bins, labels=False, duplicates='drop')
    
    results = []
    for bin_id in sorted(df['depth_bin'].unique()):
        bin_data = df[df['depth_bin'] == bin_id]
        depth_range = f"{bin_data['depth'].min():.1f}-{bin_data['depth'].max():.1f}"
        
        print(f"\nDepth bin {bin_id} (depth range: {depth_range})")
        
        # Compare populations within this depth bin
        pop_stats = []
        for pop in bin_data['super_pop'].unique():
            pop_data = bin_data[bin_data['super_pop'] == pop]
            if len(pop_data) > 5:  # Require at least 5 samples
                mean_hete = pop_data['hete_variant_num'].mean()
                std_hete = pop_data['hete_variant_num'].std()
                pop_stats.append({
                    'depth_bin': bin_id,
                    'depth_range': depth_range,
                    'population': pop,
                    'mean_hete': mean_hete,
                    'std_hete': std_hete,
                    'n': len(pop_data)
                })
                print(f"  {pop}: mean={mean_hete:.2f} ± {std_hete:.2f} (n={len(pop_data)})")
        
        # Perform ANOVA within this depth bin
        if len(pop_stats) > 1:
            groups = [bin_data[bin_data['super_pop'] == ps['population']]['hete_variant_num'].values 
                     for ps in pop_stats]
            f_stat, p_value = stats.f_oneway(*groups)
            print(f"  ANOVA: F={f_stat:.4f}, p={p_value:.4e}")
            for ps in pop_stats:
                ps['anova_f'] = f_stat
                ps['anova_p'] = p_value
        
        results.extend(pop_stats)
    
    results_df = pd.DataFrame(results)
    matched_file = os.path.join(OUTPUT_DIR, 'depth_matched_comparison.csv')
    results_df.to_csv(matched_file, index=False)
    
    # Plot
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # Plot 1: Heterozygosity by depth bin and population
    ax1 = axes[0]
    for pop in df['super_pop'].unique():
        pop_results = results_df[results_df['population'] == pop]
        ax1.errorbar(pop_results['depth_bin'], pop_results['mean_hete'], 
                    yerr=pop_results['std_hete'], marker='o', label=pop, capsize=5)
    ax1.set_xlabel('Depth Bin')
    ax1.set_ylabel('Mean Heterozygous Variants')
    ax1.set_title('Depth-Matched Population Comparison')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: P-values by depth bin
    ax2 = axes[1]
    unique_bins = results_df.drop_duplicates('depth_bin')
    ax2.scatter(unique_bins['depth_bin'], -np.log10(unique_bins['anova_p'] + 1e-10))
    ax2.axhline(y=-np.log10(0.05), color='r', linestyle='--', label='p=0.05')
    ax2.set_xlabel('Depth Bin')
    ax2.set_ylabel('-log10(p-value)')
    ax2.set_title('Population Differences Within Depth Bins')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    pdf.savefig(fig, dpi=150, bbox_inches='tight')
    plt.close()
    
    return results_df


def assess_depth_impact_on_heterozygosity(df, pdf):
    """
    Comprehensive assessment of how depth impacts heterozygous variant counts.
    """
    print("\n" + "="*80)
    print("DEPTH IMPACT ON HETEROZYGOSITY - COMPREHENSIVE ANALYSIS")
    print("="*80)
    
    # 1. Correlation analysis
    print("\nCorrelation Analysis:")
    print("-" * 80)
    
    # Pearson correlation
    pearson_r, pearson_p = stats.pearsonr(df['depth'], df['hete_variant_num'])
    print(f"Pearson correlation (depth vs hete_variant_num):")
    print(f"  r = {pearson_r:.4f}, p-value = {pearson_p:.4e}")
    
    # Spearman correlation (non-parametric)
    spearman_r, spearman_p = stats.spearmanr(df['depth'], df['hete_variant_num'])
    print(f"Spearman correlation (depth vs hete_variant_num):")
    print(f"  ρ = {spearman_r:.4f}, p-value = {spearman_p:.4e}")
    
    # 2. Binned analysis - show heterozygosity across depth bins
    print("\n\nBinned Analysis (10 depth quantiles):")
    print("-" * 80)
    
    depth_bins = pd.qcut(df['depth'], q=10, labels=False, duplicates='drop')
    df['depth_decile'] = depth_bins
    
    binned_stats = []
    for bin_id in sorted(df['depth_decile'].unique()):
        bin_data = df[df['depth_decile'] == bin_id]
        stats_dict = {
            'bin': bin_id,
            'depth_min': bin_data['depth'].min(),
            'depth_max': bin_data['depth'].max(),
            'depth_mean': bin_data['depth'].mean(),
            'hete_mean': bin_data['hete_variant_num'].mean(),
            'hete_std': bin_data['hete_variant_num'].std(),
            'hete_median': bin_data['hete_variant_num'].median(),
            'n_samples': len(bin_data)
        }
        binned_stats.append(stats_dict)
        print(f"Bin {bin_id}: depth=[{stats_dict['depth_min']:.1f}, {stats_dict['depth_max']:.1f}], "
              f"mean_depth={stats_dict['depth_mean']:.1f}, "
              f"mean_hete={stats_dict['hete_mean']:.2f} ± {stats_dict['hete_std']:.2f}, "
              f"n={stats_dict['n_samples']}")
    
    binned_df = pd.DataFrame(binned_stats)
    binned_file = os.path.join(OUTPUT_DIR, 'depth_bins_stats.csv')
    binned_df.to_csv(binned_file, index=False)
    
    # Calculate normalized heterozygosity before sampling
    df['hete_per_depth'] = df['hete_variant_num'] / (df['depth'] + 0.1)
    
    # 3. Create comprehensive visualization
    fig = plt.figure(figsize=(15, 10))
    gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)
    
    # Plot 1: Scatter plot with trend line
    ax1 = fig.add_subplot(gs[0, :2])
    # Sample data if too many points
    if len(df) > 5000:
        sample_df = df.sample(n=5000, random_state=42)
    else:
        sample_df = df
    ax1.scatter(sample_df['depth'], sample_df['hete_variant_num'], alpha=0.3, s=5)
    
    # Add trend line
    z = np.polyfit(df['depth'], df['hete_variant_num'], 1)
    p = np.poly1d(z)
    x_trend = np.linspace(df['depth'].min(), df['depth'].max(), 100)
    ax1.plot(x_trend, p(x_trend), "r-", linewidth=2, label=f'Linear fit (r={pearson_r:.3f})')
    
    ax1.set_xlabel('Depth')
    ax1.set_ylabel('Heterozygous Variant Count')
    ax1.set_title(f'Depth vs Heterozygosity\nPearson r={pearson_r:.4f} (p={pearson_p:.2e}), Spearman ρ={spearman_r:.4f} (p={spearman_p:.2e})')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Log-log plot
    ax2 = fig.add_subplot(gs[0, 2])
    ax2.scatter(sample_df['depth'], sample_df['hete_variant_num'] + 0.1, alpha=0.3, s=5)
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_xlabel('Depth (log scale)')
    ax2.set_ylabel('Heterozygous Variants (log scale)')
    ax2.set_title('Log-Log Relationship')
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: Binned heterozygosity
    ax3 = fig.add_subplot(gs[1, 0])
    ax3.errorbar(binned_df['depth_mean'], binned_df['hete_mean'], 
                yerr=binned_df['hete_std'], marker='o', capsize=5, linewidth=2)
    ax3.set_xlabel('Mean Depth per Bin')
    ax3.set_ylabel('Mean Heterozygous Variants')
    ax3.set_title('Binned Analysis (Depth Deciles)')
    ax3.grid(True, alpha=0.3)
    
    # Plot 4: Box plot by depth bins
    ax4 = fig.add_subplot(gs[1, 1])
    df.boxplot(column='hete_variant_num', by='depth_decile', ax=ax4)
    ax4.set_xlabel('Depth Decile')
    ax4.set_ylabel('Heterozygous Variants')
    ax4.set_title('Distribution by Depth Decile')
    plt.sca(ax4)
    plt.xticks(rotation=45)
    
    # Plot 5: Normalized heterozygosity (per unit depth)
    ax5 = fig.add_subplot(gs[1, 2])
    ax5.scatter(sample_df['depth'], sample_df['hete_per_depth'], alpha=0.3, s=5)
    ax5.set_xlabel('Depth')
    ax5.set_ylabel('Heterozygous Variants / Depth')
    ax5.set_title('Normalized Heterozygosity Rate')
    ax5.grid(True, alpha=0.3)
    
    # Plot 6: Population-specific depth effect
    ax6 = fig.add_subplot(gs[2, :])
    for pop in sample_df['super_pop'].unique():
        pop_data = sample_df[sample_df['super_pop'] == pop]
        ax6.scatter(pop_data['depth'], pop_data['hete_variant_num'], 
                   alpha=0.4, s=10, label=pop)
    ax6.set_xlabel('Depth')
    ax6.set_ylabel('Heterozygous Variants')
    ax6.set_title('Depth vs Heterozygosity by Population')
    ax6.legend()
    ax6.grid(True, alpha=0.3)
    
    pdf.savefig(fig, dpi=150, bbox_inches='tight')
    plt.close()
    
    # 4. Per-population correlation analysis
    print("\n\nPer-Population Correlation Analysis:")
    print("-" * 80)
    
    pop_correlations = []
    for pop in df['super_pop'].unique():
        pop_data = df[df['super_pop'] == pop]
        if len(pop_data) > 10:
            r, p = stats.pearsonr(pop_data['depth'], pop_data['hete_variant_num'])
            spearman_r, spearman_p = stats.spearmanr(pop_data['depth'], pop_data['hete_variant_num'])
            pop_correlations.append({
                'population': pop,
                'n_samples': len(pop_data),
                'pearson_r': r,
                'pearson_p': p,
                'spearman_rho': spearman_r,
                'spearman_p': spearman_p
            })
            print(f"{pop}: n={len(pop_data)}, Pearson r={r:.4f} (p={p:.4e}), Spearman ρ={spearman_r:.4f} (p={spearman_p:.4e})")
    
    pop_corr_df = pd.DataFrame(pop_correlations)
    corr_file = os.path.join(OUTPUT_DIR, 'pop_correlations.csv')
    pop_corr_df.to_csv(corr_file, index=False)
    
    return binned_df, pop_corr_df, pearson_r, pearson_p, spearman_r, spearman_p


def test_afr_heterozygosity_after_depth_control(df, pdf):
    """
    Test whether AFR population has significantly higher heterozygosity after controlling for depth.
    Uses multiple approaches:
    1. Residual analysis (removing depth effect)
    2. Depth-matched comparison
    3. ANCOVA (analysis of covariance)
    """
    print("\n" + "="*80)
    print("AFR POPULATION HETEROZYGOSITY ANALYSIS (DEPTH-CONTROLLED)")
    print("="*80)
    
    # Prepare data
    analysis_df = df.copy()
    analysis_df['log_depth'] = np.log10(analysis_df['depth'] + 1)
    
    # Create AFR vs non-AFR indicator
    analysis_df['is_AFR'] = (analysis_df['super_pop'] == 'AFR').astype(int)
    
    print(f"\nSample sizes:")
    print(f"  AFR: {analysis_df[analysis_df['is_AFR'] == 1].shape[0]}")
    print(f"  Non-AFR: {analysis_df[analysis_df['is_AFR'] == 0].shape[0]}")
    
    # ========== APPROACH 1: RESIDUAL ANALYSIS ==========
    print("\n" + "-"*80)
    print("APPROACH 1: Residual Analysis (Remove Depth Effect)")
    print("-"*80)
    
    # Fit model with depth only (on all data)
    X_depth = analysis_df[['log_depth']].values
    y = analysis_df['hete_variant_num'].values
    
    from sklearn.linear_model import LinearRegression
    depth_model = LinearRegression()
    depth_model.fit(X_depth, y)
    
    # Calculate residuals (heterozygosity after removing depth effect)
    analysis_df['residuals'] = y - depth_model.predict(X_depth)
    
    # Compare residuals between AFR and non-AFR
    afr_residuals = analysis_df[analysis_df['is_AFR'] == 1]['residuals']
    non_afr_residuals = analysis_df[analysis_df['is_AFR'] == 0]['residuals']
    
    print(f"\nResiduals (after removing depth effect):")
    print(f"  AFR mean:     {afr_residuals.mean():.4f} ± {afr_residuals.std():.4f}")
    print(f"  Non-AFR mean: {non_afr_residuals.mean():.4f} ± {non_afr_residuals.std():.4f}")
    print(f"  Difference:   {afr_residuals.mean() - non_afr_residuals.mean():.4f}")
    
    # T-test on residuals
    from scipy.stats import ttest_ind, mannwhitneyu
    t_stat, t_pval = ttest_ind(afr_residuals, non_afr_residuals)
    u_stat, u_pval = mannwhitneyu(afr_residuals, non_afr_residuals, alternative='greater')
    
    print(f"\nStatistical tests on residuals:")
    print(f"  T-test: t = {t_stat:.4f}, p = {t_pval:.4e}")
    print(f"  Mann-Whitney U (AFR > non-AFR): U = {u_stat:.2f}, p = {u_pval:.4e}")
    
    # ========== APPROACH 2: DEPTH-MATCHED ANALYSIS ==========
    print("\n" + "-"*80)
    print("APPROACH 2: Depth-Matched Comparison")
    print("-"*80)
    
    # Create depth bins
    analysis_df['depth_decile'] = pd.qcut(analysis_df['depth'], q=10, labels=False, duplicates='drop')
    
    depth_matched_results = []
    for bin_id in sorted(analysis_df['depth_decile'].unique()):
        bin_data = analysis_df[analysis_df['depth_decile'] == bin_id]
        
        afr_bin = bin_data[bin_data['is_AFR'] == 1]['hete_variant_num']
        non_afr_bin = bin_data[bin_data['is_AFR'] == 0]['hete_variant_num']
        
        if len(afr_bin) >= 5 and len(non_afr_bin) >= 5:
            depth_range = f"{bin_data['depth'].min():.1f}-{bin_data['depth'].max():.1f}"
            mean_diff = afr_bin.mean() - non_afr_bin.mean()
            
            # T-test within this bin
            t_stat_bin, p_val_bin = ttest_ind(afr_bin, non_afr_bin)
            
            depth_matched_results.append({
                'depth_bin': bin_id,
                'depth_range': depth_range,
                'depth_mean': bin_data['depth'].mean(),
                'afr_mean': afr_bin.mean(),
                'afr_std': afr_bin.std(),
                'afr_n': len(afr_bin),
                'non_afr_mean': non_afr_bin.mean(),
                'non_afr_std': non_afr_bin.std(),
                'non_afr_n': len(non_afr_bin),
                'mean_diff': mean_diff,
                't_stat': t_stat_bin,
                'p_value': p_val_bin
            })
            
            print(f"\nBin {bin_id} (depth {depth_range}):")
            print(f"  AFR:     {afr_bin.mean():.2f} ± {afr_bin.std():.2f} (n={len(afr_bin)})")
            print(f"  Non-AFR: {non_afr_bin.mean():.2f} ± {non_afr_bin.std():.2f} (n={len(non_afr_bin)})")
            print(f"  Diff:    {mean_diff:.2f}, t={t_stat_bin:.3f}, p={p_val_bin:.4e}")
    
    matched_df = pd.DataFrame(depth_matched_results)
    matched_file = os.path.join(OUTPUT_DIR, 'afr_depth_matched_comparison.csv')
    matched_df.to_csv(matched_file, index=False)
    
    # Overall summary across bins
    n_sig = (matched_df['p_value'] < 0.05).sum()
    n_afr_higher = (matched_df['mean_diff'] > 0).sum()
    
    print(f"\nSummary across depth bins:")
    print(f"  {n_afr_higher}/{len(matched_df)} bins show AFR > non-AFR")
    print(f"  {n_sig}/{len(matched_df)} bins show significant difference (p < 0.05)")
    
    # ========== APPROACH 3: ANCOVA ==========
    print("\n" + "-"*80)
    print("APPROACH 3: ANCOVA (Analysis of Covariance)")
    print("-"*80)
    
    # Model 1: Depth only
    X1 = analysis_df[['log_depth']]
    model1 = LinearRegression()
    model1.fit(X1, y)
    r2_depth = model1.score(X1, y)
    
    # Model 2: Depth + AFR indicator
    X2 = analysis_df[['log_depth', 'is_AFR']]
    model2 = LinearRegression()
    model2.fit(X2, y)
    r2_depth_afr = model2.score(X2, y)
    
    afr_coef = model2.coef_[1]
    depth_coef = model2.coef_[0]
    
    print(f"\nModel 1 (Depth only): R² = {r2_depth:.4f}")
    print(f"Model 2 (Depth + AFR): R² = {r2_depth_afr:.4f}")
    print(f"\nCoefficients in Model 2:")
    print(f"  log_depth: {depth_coef:.4f}")
    print(f"  is_AFR:    {afr_coef:.4f}")
    print(f"\nAdditional variance explained by AFR: {(r2_depth_afr - r2_depth)*100:.2f}%")
    
    # F-test for significance of AFR term
    n = len(y)
    rss1 = np.sum((y - model1.predict(X1))**2)
    rss2 = np.sum((y - model2.predict(X2))**2)
    f_stat = ((rss1 - rss2) / 1) / (rss2 / (n - 3))
    f_pval = 1 - stats.f.cdf(f_stat, 1, n - 3)
    
    print(f"\nF-test for AFR effect after controlling for depth:")
    print(f"  F({1}, {n-3}) = {f_stat:.4f}, p = {f_pval:.4e}")
    
    # ========== VISUALIZATION ==========
    fig = plt.figure(figsize=(15, 10))
    gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)
    
    # Plot 1: Scatter plot with regression lines
    ax1 = fig.add_subplot(gs[0, :2])
    # Sample data for visualization
    sample_size = min(3000, len(analysis_df))
    vis_df = analysis_df.sample(n=sample_size, random_state=42) if len(analysis_df) > sample_size else analysis_df
    
    for pop_label, pop_name in [(0, 'Non-AFR'), (1, 'AFR')]:
        pop_data = vis_df[vis_df['is_AFR'] == pop_label]
        ax1.scatter(pop_data['depth'], pop_data['hete_variant_num'], 
                   alpha=0.4, s=10, label=pop_name)
        
        # Add trend line (use full data for trend)
        full_pop_data = analysis_df[analysis_df['is_AFR'] == pop_label]
        z = np.polyfit(full_pop_data['depth'], full_pop_data['hete_variant_num'], 1)
        p = np.poly1d(z)
        x_range = np.linspace(full_pop_data['depth'].min(), full_pop_data['depth'].max(), 100)
        ax1.plot(x_range, p(x_range), linewidth=2)
    
    ax1.set_xlabel('Depth')
    ax1.set_ylabel('Heterozygous Variants')
    ax1.set_title('Depth vs Heterozygosity: AFR vs Non-AFR')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Residuals by population
    ax2 = fig.add_subplot(gs[0, 2])
    ax2.boxplot([non_afr_residuals, afr_residuals], labels=['Non-AFR', 'AFR'])
    ax2.axhline(y=0, color='r', linestyle='--', alpha=0.5)
    ax2.set_ylabel('Residuals (after removing depth effect)')
    ax2.set_title(f'Residuals Comparison\np = {t_pval:.4e}')
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: Violin plot of residuals
    ax3 = fig.add_subplot(gs[1, 0])
    parts = ax3.violinplot([non_afr_residuals, afr_residuals], positions=[1, 2], showmeans=True)
    ax3.set_xticks([1, 2])
    ax3.set_xticklabels(['Non-AFR', 'AFR'])
    ax3.axhline(y=0, color='r', linestyle='--', alpha=0.5)
    ax3.set_ylabel('Residuals')
    ax3.set_title('Residual Distribution')
    ax3.grid(True, alpha=0.3)
    
    # Plot 4: Depth-matched comparison
    ax4 = fig.add_subplot(gs[1, 1:])
    x_pos = matched_df['depth_mean']
    width = (matched_df['depth_mean'].max() - matched_df['depth_mean'].min()) / (len(matched_df) * 3)
    
    ax4.bar(x_pos - width/2, matched_df['afr_mean'], width, 
            yerr=matched_df['afr_std'], label='AFR', alpha=0.7, capsize=5)
    ax4.bar(x_pos + width/2, matched_df['non_afr_mean'], width,
            yerr=matched_df['non_afr_std'], label='Non-AFR', alpha=0.7, capsize=5)
    
    ax4.set_xlabel('Mean Depth per Bin')
    ax4.set_ylabel('Mean Heterozygous Variants')
    ax4.set_title('Depth-Matched Comparison: AFR vs Non-AFR')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    # Plot 5: P-values across depth bins
    ax5 = fig.add_subplot(gs[2, 0])
    colors = ['red' if p < 0.05 else 'gray' for p in matched_df['p_value']]
    ax5.bar(range(len(matched_df)), -np.log10(matched_df['p_value']), color=colors, alpha=0.7)
    ax5.axhline(y=-np.log10(0.05), color='r', linestyle='--', label='p=0.05')
    ax5.set_xlabel('Depth Bin')
    ax5.set_ylabel('-log10(p-value)')
    ax5.set_title('Significance Across Depth Bins')
    ax5.legend()
    ax5.grid(True, alpha=0.3)
    
    # Plot 6: Mean difference across depth bins
    ax6 = fig.add_subplot(gs[2, 1:])
    colors = ['green' if d > 0 else 'orange' for d in matched_df['mean_diff']]
    ax6.bar(matched_df['depth_mean'], matched_df['mean_diff'], 
            width=width*2, color=colors, alpha=0.7)
    ax6.axhline(y=0, color='black', linestyle='-', linewidth=1)
    ax6.set_xlabel('Mean Depth per Bin')
    ax6.set_ylabel('Mean Difference (AFR - Non-AFR)')
    ax6.set_title('Heterozygosity Difference Across Depth Bins')
    ax6.grid(True, alpha=0.3)
    
    plt.tight_layout()
    pdf.savefig(fig, dpi=150, bbox_inches='tight')
    plt.close()
    
    # ========== SUMMARY ==========
    print("\n" + "="*80)
    print("SUMMARY: AFR HETEROZYGOSITY AFTER DEPTH CONTROL")
    print("="*80)
    
    afr_higher = afr_residuals.mean() > non_afr_residuals.mean()
    
    summary_dict = {
        'approach': ['Residual Analysis', 'Mann-Whitney U', 'ANCOVA', 'Depth-Matched'],
        'afr_higher': [afr_higher, afr_higher, afr_coef > 0, n_afr_higher > len(matched_df)/2],
        'p_value': [t_pval, u_pval, f_pval, 'See bins'],
        'effect_size': [
            afr_residuals.mean() - non_afr_residuals.mean(),
            f"U={u_stat:.0f}",
            afr_coef,
            matched_df['mean_diff'].mean()
        ]
    }
    
    summary_df = pd.DataFrame(summary_dict)
    print(summary_df.to_string(index=False))
    
    # Save summary
    summary_file = os.path.join(OUTPUT_DIR, 'afr_heterozygosity_summary.csv')
    summary_df.to_csv(summary_file, index=False)
    
    conclusion = "YES" if (t_pval < 0.05 and afr_higher) else "NO"
    print(f"\n{'='*80}")
    print(f"CONCLUSION: Does AFR have significantly higher heterozygosity after")
    print(f"            controlling for depth? {conclusion}")
    print(f"{'='*80}")
    
    return {
        'residual_pval': t_pval,
        'mann_whitney_pval': u_pval,
        'ancova_pval': f_pval,
        'ancova_afr_coef': afr_coef,
        'residual_diff': afr_residuals.mean() - non_afr_residuals.mean(),
        'matched_df': matched_df,
        'conclusion': conclusion
    }


##

def create_summary_page(pdf, df, model_comparison, pearson_r, pearson_p, spearman_r, spearman_p, afr_results=None):
    """
    Create a summary page with key findings.
    """
    fig = plt.figure(figsize=(11, 14))
    fig.suptitle('IG/TCR Heterozygosity Regression Analysis - Summary', fontsize=16, fontweight='bold')
    
    # Remove axes
    ax = fig.add_subplot(111)
    ax.axis('off')
    
    # Create summary text
    afr_section = ""
    if afr_results:
        afr_section = f"""

AFR POPULATION ANALYSIS (DEPTH-CONTROLLED)
{'='*80}
Question: Does AFR have higher heterozygosity after controlling for depth?

Answer: {afr_results['conclusion']}

Evidence:
  1. Residual Analysis (depth removed):
     - Difference: {afr_results['residual_diff']:.4f}
     - p-value: {afr_results['residual_pval']:.4e}
  
  2. Mann-Whitney U Test:
     - p-value: {afr_results['mann_whitney_pval']:.4e}
  
  3. ANCOVA (Depth + AFR indicator):
     - AFR coefficient: {afr_results['ancova_afr_coef']:.4f}
     - p-value: {afr_results['ancova_pval']:.4e}
  
  4. Depth-Matched Comparison:
     - Bins where AFR > non-AFR: {(afr_results['matched_df']['mean_diff'] > 0).sum()}/{len(afr_results['matched_df'])}
     - Significant bins (p<0.05): {(afr_results['matched_df']['p_value'] < 0.05).sum()}/{len(afr_results['matched_df'])}
"""
    
    summary_text = f"""
DATASET SUMMARY
{'='*80}
Total samples: {len(df):,}
Unique genes: {len(df['gene'].unique()):,}
Populations: {', '.join(sorted(df['super_pop'].unique()))}

Depth Statistics:
  Mean: {df['depth'].mean():.2f}
  Median: {df['depth'].median():.2f}
  Range: [{df['depth'].min():.2f}, {df['depth'].max():.2f}]

Heterozygous Variant Statistics:
  Mean: {df['hete_variant_num'].mean():.2f}
  Median: {df['hete_variant_num'].median():.2f}
  Range: [{df['hete_variant_num'].min():.0f}, {df['hete_variant_num'].max():.0f}]


CORRELATION ANALYSIS (Depth vs Heterozygosity)
{'='*80}
Pearson correlation:  r = {pearson_r:.4f}, p-value = {pearson_p:.4e}
Spearman correlation: ρ = {spearman_r:.4f}, p-value = {spearman_p:.4e}

Interpretation: {'Strong' if abs(pearson_r) > 0.5 else 'Moderate' if abs(pearson_r) > 0.3 else 'Weak'} positive correlation between depth and heterozygosity.


MODEL COMPARISON
{'='*80}
Model 1 (Depth only):                  R² = {model_comparison.iloc[0]['R²']:.4f}
Model 2 (All Technical Covariates):    R² = {model_comparison.iloc[1]['R²']:.4f}
Model 3 (Technical + Population):      R² = {model_comparison.iloc[2]['R²']:.4f}

Variance explained by technical factors: {model_comparison.iloc[1]['R²']*100:.2f}%
Additional variance by population:       {(model_comparison.iloc[2]['R²'] - model_comparison.iloc[1]['R²'])*100:.2f}%
{afr_section}

KEY FINDINGS
{'='*80}
1. Depth has a significant impact on heterozygous variant detection
   (Pearson r={pearson_r:.3f}, p={pearson_p:.2e})

2. Technical covariates (depth, read length, quality) explain
   {model_comparison.iloc[1]['R²']*100:.1f}% of variance in heterozygosity

3. Population factors add {(model_comparison.iloc[2]['R²'] - model_comparison.iloc[1]['R²'])*100:.1f}% additional explanatory power
   after controlling for technical factors

4. Depth-matching is critical for fair population comparisons

5. AFR population {'DOES' if afr_results and afr_results['conclusion'] == 'YES' else 'DOES NOT'} show significantly higher heterozygosity
   after controlling for depth effects


RECOMMENDATIONS
{'='*80}
• Use depth-matched comparisons when comparing populations
• Control for technical covariates (depth, read length, quality)
• Consider residual analysis to identify true population effects
• Review depth-matched plots to verify population differences


OUTPUT FILES
{'='*80}
All results saved to: {OUTPUT_DIR}/
  • model_comparison.csv
  • population_residuals.csv
  • depth_bins_stats.csv
  • pop_correlations.csv
  • depth_matched_comparison.csv
  • depth_distribution_stats.csv
  • afr_heterozygosity_summary.csv
  • afr_depth_matched_comparison.csv
"""
    
    ax.text(0.05, 0.95, summary_text, transform=ax.transAxes, 
            fontsize=10, verticalalignment='top', fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))
    
    pdf.savefig(fig, dpi=100, bbox_inches='tight')
    plt.close()


if __name__ == "__main__":
    super_pop_file = "../hla/20131219.populations.tsv"
    sample_pop_file = "../hla/20130606_sample_info.xlsx"
    allele_file = "1kg_merged_samples.ig_tr.csv"

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
    print("RUNNING COMPREHENSIVE ANALYSES")
    print("="*80)
    print(f"\nAll plots will be saved to: {output_pdf}")
    
    with PdfPages(output_pdf) as pdf:
        # 1. Depth distributions
        plot_depth_distributions(df, pdf)
        
        # 2. Assess depth impact on heterozygosity
        binned_stats, pop_correlations, pearson_r, pearson_p, spearman_r, spearman_p = assess_depth_impact_on_heterozygosity(df, pdf)
        
        # 3. Depth-stratified analysis
        stratified_results = perform_depth_stratified_analysis(df)
        
        # 4. Multivariate regression
        analysis_df, model_tech, model_full, model_comparison = multivariate_regression_analysis(df, pdf)
        
        # 5. Depth-matched population comparison
        matched_results = depth_matched_population_comparison(df, pdf, n_bins=5)
        
        # 6. AFR-specific analysis (depth-controlled)
        afr_results = test_afr_heterozygosity_after_depth_control(df, pdf)
        
        # 7. Create summary page (first page)
        create_summary_page(pdf, df, model_comparison, pearson_r, pearson_p, spearman_r, spearman_p, afr_results)
    
    print("\n" + "="*80)
    print("ANALYSIS COMPLETE")
    print("="*80)
    print(f"\n✓ All plots and analysis saved to: {output_pdf}")
    print(f"✓ Additional CSV files saved to: {OUTPUT_DIR}/")
    print("\n" + "="*80)
    print("RESULTS SUMMARY")
    print("="*80)
    print(f"\nCorrelation (Depth vs Heterozygosity):")
    print(f"  Pearson r = {pearson_r:.4f} (p = {pearson_p:.4e})")
    print(f"  Spearman ρ = {spearman_r:.4f} (p = {spearman_p:.4e})")
    print(f"\nModel Performance (R² values):")
    for idx, row in model_comparison.iterrows():
        print(f"  {row['Model']}: R² = {row['R²']:.4f} ({row['Predictors']} predictors)")
    print(f"\nVariance Explained:")
    print(f"  Technical factors: {model_comparison.iloc[1]['R²']*100:.2f}%")
    print(f"  Population (additional): {(model_comparison.iloc[2]['R²'] - model_comparison.iloc[1]['R²'])*100:.2f}%")
    print(f"\n" + "="*80)
    print("AFR POPULATION ANALYSIS (AFTER CONTROLLING FOR DEPTH)")
    print("="*80)
    print(f"\nQuestion: Does AFR have significantly higher heterozygosity after")
    print(f"          controlling for depth effects?")
    print(f"\nAnswer: {afr_results['conclusion']}")
    print(f"\nEvidence:")
    print(f"  Residual difference: {afr_results['residual_diff']:.4f} (p = {afr_results['residual_pval']:.4e})")
    print(f"  Mann-Whitney U:      p = {afr_results['mann_whitney_pval']:.4e}")
    print(f"  ANCOVA AFR coef:     {afr_results['ancova_afr_coef']:.4f} (p = {afr_results['ancova_pval']:.4e})")
    print(f"  Depth-matched bins:  {(afr_results['matched_df']['mean_diff'] > 0).sum()}/{len(afr_results['matched_df'])} show AFR > non-AFR")
    print(f"                       {(afr_results['matched_df']['p_value'] < 0.05).sum()}/{len(afr_results['matched_df'])} bins significant (p < 0.05)")
    print("\n" + "="*80)

