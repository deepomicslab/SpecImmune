import os
import subprocess
import numpy as np
import argparse
from scipy.spatial.distance import squareform
from tqdm import tqdm
import pandas as pd
import random

# 基因列表
GENES = [
    'HLA-A', 'HLA-B', 'HLA-C',
    'HLA-DPA1', 'HLA-DPB1', 'HLA-DPB2', 'HLA-DQA1', 'HLA-DQB1',
    'HLA-DRB1', 'HLA-DQA2'
]


def generate_bam_and_get_divergency(allele1, allele2, bam_file):
    """使用 minimap2 比对两个等位基因，并返回比对的差异值。"""
    minimap2_cmd = f"minimap2 -ax asm10 {allele1} {allele2} -t {args.threads} | samtools view -bS -F 0x800 - | samtools sort > {bam_file}"
    subprocess.run(minimap2_cmd, shell=True, check=True, stderr=subprocess.DEVNULL)
    
    result = subprocess.run(f"samtools view {bam_file} | head -n 1", 
                            stdout=subprocess.PIPE, shell=True, check=True, stderr=subprocess.DEVNULL)
    
    first_record = result.stdout.decode('utf-8').strip()
    
    if not first_record:
        return None
    
    for field in first_record.split('\t'):
        if field.startswith('de:f:'):
            return float(field.split(':')[-1])
    
    return None

def generate_bam_and_get_divergency_nm(allele1, allele2, bam_file, threads=4):
    """
    使用 minimap2 比对两个等位基因，并返回比对的差异值 (使用 NM 标签)。
    
    Parameters:
    - allele1: 第一个等位基因的路径
    - allele2: 第二个等位基因的路径
    - bam_file: 输出 BAM 文件的路径
    - threads: 使用的线程数 (默认值: 4)
    
    Returns:
    - NM tag (edit distance) from the first alignment record, or None if not found.
    """
    # 运行 minimap2 进行比对并生成 BAM 文件
    minimap2_cmd = f"minimap2 -ax asm10 {allele1} {allele2} -t {threads} | samtools view -bS -F 0x800 - | samtools sort > {bam_file}"
    subprocess.run(minimap2_cmd, shell=True, check=True, stderr=subprocess.DEVNULL)
    
    # 提取 BAM 文件中的第一个比对记录
    result = subprocess.run(f"samtools view {bam_file} | head -n 1", 
                            stdout=subprocess.PIPE, shell=True, check=True, stderr=subprocess.DEVNULL)
    
    # 解析第一个比对记录
    first_record = result.stdout.decode('utf-8').strip()
    
    if not first_record:
        return None
    
    # 查找 NM 标签 (edit distance)
    for field in first_record.split('\t'):
        if field.startswith('NM:i:'):
            return int(field.split(':')[-1])  # 返回 NM 标签的值
    
    return None

def compute_sample_pair_divergency_vector(sample1, sample2, bam_output_folder):
    """计算两个样本之间的差异向量（针对每个基因），返回差异值。"""
    divergency_vector = []
    sample1_folder = os.path.join(args.path, sample1, sample1, 'Sequences')
    sample2_folder = os.path.join(args.path, sample2, sample2, 'Sequences')
    
    for gene in GENES:
        s1_allele1 = os.path.join(sample1_folder, f'HLA.allele.1.{gene}.fasta')
        s2_allele1 = os.path.join(sample2_folder, f'HLA.allele.1.{gene}.fasta')
        s1_allele2 = os.path.join(sample1_folder, f'HLA.allele.2.{gene}.fasta')
        s2_allele2 = os.path.join(sample2_folder, f'HLA.allele.2.{gene}.fasta')

        dv_s1_s2 = []
        # 比对等位基因1
        dv_s1a1_s2a1 = generate_bam_and_get_divergency_nm(s1_allele1, s2_allele1, os.path.join(bam_output_folder, f'{sample1}_vs_{sample2}_{gene}_s1a1_s2a1.bam')) or 0
        dv_s1a2_s2a2 = generate_bam_and_get_divergency_nm(s1_allele2, s2_allele2, os.path.join(bam_output_folder, f'{sample1}_vs_{sample2}_{gene}_s1a2_s2a2.bam')) or 0
        dv_s1_s2.append((dv_s1a1_s2a1 + dv_s1a2_s2a2) / 2)

        # 跨等位基因比对
        dv_s1a1_s2a2 = generate_bam_and_get_divergency_nm(s1_allele1, s2_allele2, os.path.join(bam_output_folder, f'{sample1}_vs_{sample2}_{gene}_s1a1_s2a2.bam')) or 0
        dv_s1a2_s2a1 = generate_bam_and_get_divergency_nm(s1_allele2, s2_allele1, os.path.join(bam_output_folder, f'{sample1}_vs_{sample2}_{gene}_s1a2_s2a1.bam')) or 0
        dv_s1_s2.append((dv_s1a1_s2a2 + dv_s1a2_s2a1) / 2)

        # 选择最小差异值作为该基因的差异值
        divergency_vector.append(min(dv_s1_s2))
    
    return divergency_vector

def compute_distance_matrix(sample_names, output_csv, bam_output_folder):
    """计算所有样本两两之间的距离矩阵并保存到 CSV。"""
    num_samples = len(sample_names)
    num_genes = len(GENES)
    distance_matrix = np.zeros((num_samples, num_samples, num_genes))
    
    os.makedirs(bam_output_folder, exist_ok=True)
    
    total_pairs = num_samples * (num_samples - 1) // 2
    for i in range(num_samples):
        for j in range(i + 1, num_samples):
            tqdm.write(f"Processing {sample_names[i]} vs {sample_names[j]}... ({i * num_samples + j + 1}/{total_pairs})")
            divergency_vector = compute_sample_pair_divergency_vector(sample_names[i], sample_names[j], bam_output_folder)
            distance_matrix[i, j, :] = divergency_vector
            distance_matrix[j, i, :] = divergency_vector
    
    # 保存结果到 CSV 文件
    with open(output_csv, 'w') as f:
        header = "Sample1,Sample2," + ",".join(GENES) + "\n"
        f.write(header)
        
        for i in range(num_samples):
            for j in range(i + 1, num_samples):
                sample1 = sample_names[i]
                sample2 = sample_names[j]
                divergency_vector = distance_matrix[i, j, :]
                line = f"{sample1},{sample2}," + ",".join(map(str, divergency_vector)) + "\n"
                f.write(line)

    print(f"Distance matrix saved to {output_csv}")

def get_population(df, sample_name):
    """根据样本名称获取其对应的 Population。"""
    sample_row = df[df['Sample'] == sample_name]
    
    if sample_row.empty:
        return "Unknown"
    
    return sample_row['Population'].values[0]

def parse_arguments():
    """解析命令行参数。"""
    parser = argparse.ArgumentParser(description="计算样本之间的HLA基因等位基因距离矩阵。")
    
    parser.add_argument('-i', '--input', required=True, help="包含样本名称的输入文件（每行一个样本）。")
    parser.add_argument('-p', '--path', required=True, help="样本文件夹的基本路径。")
    parser.add_argument('-o', '--output', required=True, help="输出距离矩阵的CSV文件。")
    parser.add_argument('-b', '--bam_output', required=True, help="保存生成BAM文件的目录。")
    parser.add_argument('-t', '--threads', type=int, default=1, help="minimap2使用的线程数。")
    parser.add_argument('-d', '--depth', required=True, help="包含基因深度信息的CSV文件。")
    parser.add_argument('-min_depth', '--min_depth', type=int, default=10, help="基因深度过滤的最小值。")
    parser.add_argument('-cs', '--chosen_samples_output', required=True, help="保存选择样本的CSV文件。")
    
    return parser.parse_args()

def filter_samples_by_depth(depth_df, min_depth, genes):
    """根据GENES列表中的基因深度过滤样本。"""
    # 确保第一列是 "Sample"
    if depth_df.columns[0] != 'Sample':
        depth_df.rename(columns={depth_df.columns[0]: 'Sample'}, inplace=True)

    # 只选择GENES列表中的基因列进行过滤
    selected_columns = ['Sample'] + [gene for gene in genes if gene in depth_df.columns]
    
    # 过滤掉深度不足的样本
    filtered_df = depth_df[selected_columns]
    filtered_df = filtered_df[(filtered_df.iloc[:, 1:] >= min_depth).all(axis=1)]
    
    return filtered_df['Sample'].unique()

def choose_random_samples_by_superpopulation(sample_df, num_samples_per_group=10):
    """
    根据Superpopulation随机选择样本，最多选择 num_samples_per_group 个样本。
    如果Superpopulation的样本数少于 num_samples_per_group，则选择所有样本。
    """
    grouped = sample_df.groupby('Superpopulation')
    
    selected_samples = []
    
    for superpop, group in grouped:
        # 从每个 Superpopulation 中选择最多 num_samples_per_group 个样本
        selected_samples_from_group = group.sample(n=min(len(group), num_samples_per_group))['Sample'].values.tolist()
        selected_samples.extend(selected_samples_from_group)
    
    return selected_samples

if __name__ == '__main__':
    args = parse_arguments()
    
    # 读取输入样本文件
    sample_df = pd.read_csv(args.input, header=None, names=["Sample"])
    metadata_df = pd.read_excel('./20130606_sample_info.xlsx', engine='openpyxl')

    # 处理列名，去掉空格
    metadata_df.columns = metadata_df.columns.str.replace(' ', '_')

    # 获取样本对应的Population
    sample_df['Population'] = sample_df['Sample'].apply(lambda x: get_population(metadata_df, x))
    
    # Population到Superpopulation的映射字典
    population_superpopulation_dict = {
        'CDX': 'EAS', 'CHB': 'EAS', 'JPT': 'EAS', 'KHV': 'EAS', 'CHS': 'EAS',
        'BEB': 'SAS', 'GIH': 'SAS', 'ITU': 'SAS', 'PJL': 'SAS', 'STU': 'SAS',
        'ASW': 'AFR', 'ACB': 'AFR', 'ESN': 'AFR', 'GWD': 'AFR', 'LWK': 'AFR',
        'MSL': 'AFR', 'YRI': 'AFR', 'GBR': 'EUR', 'FIN': 'EUR', 'IBS': 'EUR',
        'TSI': 'EUR', 'CEU': 'EUR', 'CLM': 'AMR', 'MXL': 'AMR', 'PEL': 'AMR', 
        'PUR': 'AMR'
    }

    # 映射Population到Superpopulation
    sample_df['Superpopulation'] = sample_df['Population'].map(population_superpopulation_dict)

    # 读取深度文件
    depth_df = pd.read_csv(args.depth)

    # 检查并重命名第一列为 'Sample' 如果它不是 'Sample'
    if depth_df.columns[0] != 'Sample':
        depth_df.rename(columns={depth_df.columns[0]: 'Sample'}, inplace=True)

    # 根据GENES列表中的基因深度过滤样本
    filtered_samples = filter_samples_by_depth(depth_df, args.min_depth, GENES)

    # 过滤后的样本
    filtered_sample_df = sample_df[sample_df['Sample'].isin(filtered_samples)]

    # 随机选择每个Superpopulation的样本（尽可能选择10个样本）
    selected_samples = choose_random_samples_by_superpopulation(filtered_sample_df, num_samples_per_group=10)

    # 从DataFrame中获取选择的样本
    selected_samples_df = sample_df[sample_df['Sample'].isin(selected_samples)]

    # 保存选择的样本到CSV文件
    selected_samples_df.to_csv(args.chosen_samples_output, index=False)
    print(f"Selected samples saved to {args.chosen_samples_output}")

    # 计算选择样本的距离矩阵
    compute_distance_matrix(selected_samples, args.output, args.bam_output)