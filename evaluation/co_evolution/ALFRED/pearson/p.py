import pandas as pd
from scipy.stats import pearsonr
import sys

def read_gene_frequency_file(filename):
    """
    读取基因频率文件并返回DataFrame
    """
    try:
        # 读取文件，跳过分隔线
        with open(filename, 'r') as f:
            lines = f.readlines()
        
        # 找到数据开始的行（跳过标题和分隔线）
        data_start = 0
        for i, line in enumerate(lines):
            if line.strip().startswith('Population|'):
                data_start = i + 2  # 跳过标题行和分隔线
                break
        
        # 提取数据行
        data_lines = []
        for line in lines[data_start:]:
            line = line.strip()
            if line and not line.startswith('-'):
                data_lines.append(line)
        
        # 解析数据
        data = []
        for line in data_lines:
            parts = line.split('|')
            if len(parts) >= 6:
                population = parts[0].strip()
                sample_id = parts[1].strip()
                date = parts[2].strip()
                n2 = parts[3].strip()
                ref_freq = float(parts[4].strip())
                alt_freq = float(parts[5].strip())
                
                data.append({
                    'Population': population,
                    'SampleUId': sample_id,
                    'Date': date,
                    '2N': n2,
                    'ref_freq': ref_freq,
                    'alt_freq': alt_freq
                })
        
        return pd.DataFrame(data)
    
    except Exception as e:
        print(f"读取文件 {filename} 时出错: {e}")
        return None

def calculate_correlation(file1, file2, use_ref=True):
    """
    计算两个基因位点的皮尔森相关系数
    
    Parameters:
    file1: 第一个基因位点文件路径
    file2: 第二个基因位点文件路径
    use_ref: True使用ref频率，False使用alt频率
    """
    
    # 读取两个文件
    df1 = read_gene_frequency_file(file1)
    df2 = read_gene_frequency_file(file2)
    
    if df1 is None or df2 is None:
        print("文件读取失败")
        return None, None
    
    print(f"文件1包含 {len(df1)} 个样本")
    print(f"文件2包含 {len(df2)} 个样本")
    
    # 基于Population和SampleUId匹配样本
    merged = pd.merge(df1, df2, on=['Population', 'SampleUId'], suffixes=('_locus1', '_locus2'))
    
    if len(merged) == 0:
        print("没有找到匹配的样本")
        return None, None
    
    print(f"找到 {len(merged)} 个匹配的样本")
    
    # 选择要比较的频率列
    freq_col = 'ref_freq' if use_ref else 'alt_freq'
    locus1_freqs = merged[f'{freq_col}_locus1'].values
    locus2_freqs = merged[f'{freq_col}_locus2'].values
    
    # 计算皮尔森相关系数
    corr, p_value = pearsonr(locus1_freqs, locus2_freqs)
    
    return corr, p_value, merged

def main():
    """
    主函数
    """
    if len(sys.argv) != 3:
        print("使用方法: python script.py <文件1路径> <文件2路径>")
        print("例如: python script.py locus1.txt locus2.txt")
        return
    
    file1 = sys.argv[1]
    file2 = sys.argv[2]
    
    print("="*60)
    print("基因位点频率相关性分析")
    print("="*60)
    
    # 计算ref频率的相关性
    print("\n分析ref碱基频率相关性...")
    corr_ref, p_ref, merged_data = calculate_correlation(file1, file2, use_ref=True)
    
    if corr_ref is not None:
        print(f"ref频率皮尔森相关系数: {corr_ref:.4f}")
        print(f"p值: {p_ref:.4e}")
        print(f"相关性强度: ", end="")
        if abs(corr_ref) >= 0.7:
            print("强相关")
        elif abs(corr_ref) >= 0.3:
            print("中等相关")
        else:
            print("弱相关")
    
    # 计算alt频率的相关性
    print("\n分析alt碱基频率相关性...")
    corr_alt, p_alt, _ = calculate_correlation(file1, file2, use_ref=False)
    
    if corr_alt is not None:
        print(f"alt频率皮尔森相关系数: {corr_alt:.4f}")
        print(f"p值: {p_alt:.4e}")
        print(f"相关性强度: ", end="")
        if abs(corr_alt) >= 0.7:
            print("强相关")
        elif abs(corr_alt) >= 0.3:
            print("中等相关")
        else:
            print("弱相关")
    
    # 显示匹配的样本信息
    if merged_data is not None and len(merged_data) > 0:
        print(f"\n匹配的群体分布:")
        population_counts = merged_data['Population'].value_counts()
        for pop, count in population_counts.items():
            print(f"  {pop}: {count} 个样本")

if __name__ == "__main__":
    main()