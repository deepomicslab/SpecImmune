import os
import subprocess
import numpy as np
import argparse
from scipy.spatial.distance import squareform
from tqdm import tqdm

GENES = [ 'HLA-A', 'HLA-B', 'HLA-C', 'HLA-DMA', 'HLA-DMB', 'HLA-DOA', 'HLA-DOB', 'HLA-DPA1', 'HLA-DPB1', 'HLA-DPB2', 'HLA-DQA1', 'HLA-DQB1','HLA-DQB2', 'HLA-DRA', 'HLA-DRB1', 'HLA-DRB3', 'HLA-DRB4', 'HLA-DRB5', 'HLA-E', 'HLA-F', 'HLA-G', 'HLA-H', 'HLA-J', 'HLA-K', 'HLA-L', 'HLA-P', 'HLA-V', 'HLA-Y', 'HLA-DQA2', 'HLA-DPA2', 'HLA-N', 'HLA-S', 'HLA-T', 'HLA-U', 'HLA-W', 'MICA', 'MICB', 'TAP1', 'TAP2', 'HFE' ]

def generate_bam_and_get_divergency(allele1, allele2, bam_file):
    # 使用 minimap2 和 samtools 生成 BAM 文件，同时抑制 samtools 输出日志信息
    minimap2_cmd = f"minimap2 -ax asm10 {allele1} {allele2} -t {args.threads} | samtools view -bS -F 0x800 - | samtools sort > {bam_file}"
    subprocess.run(minimap2_cmd, shell=True, check=True, stderr=subprocess.DEVNULL)
    
    # 使用 samtools 提取 BAM 文件中的第一条记录，同时抑制 samtools 输出日志信息
    result = subprocess.run(f"samtools view {bam_file} | head -n 1", 
                            stdout=subprocess.PIPE, shell=True, check=True, stderr=subprocess.DEVNULL)
    
    first_record = result.stdout.decode('utf-8').strip()
    
    if not first_record:
        return None  # 如果没有任何记录，返回 None
    
    # 查找 de:f: 标签并提取 divergency 值
    for field in first_record.split('\t'):
        if field.startswith('de:f:'):
            return float(field.split(':')[-1])
    
    # 如果没有找到 de:f: 标签，返回 None
    return None

def compute_sample_pair_divergency_vector(sample1, sample2, bam_output_folder):
    divergency_vector = []
    sample1_folder = os.path.join(args.path, sample1, sample1, 'Sequences')
    sample2_folder = os.path.join(args.path, sample2, sample2, 'Sequences')
    
    for gene in GENES:
        s1_allele1 = os.path.join(sample1_folder, f'HLA.allele.1.{gene}.fasta')
        s2_allele1 = os.path.join(sample2_folder, f'HLA.allele.1.{gene}.fasta')
        s1_allele2 = os.path.join(sample1_folder, f'HLA.allele.2.{gene}.fasta')
        s2_allele2 = os.path.join(sample2_folder, f'HLA.allele.2.{gene}.fasta')

        dv_s1_s2 = []
        dv_s1a1_s2a1 = 0
        dv_s1a2_s2a2 = 0
        dv_s1a1_s2a2 = 0
        dv_s1a2_s2a1 = 0
        # for s1a1-s2a1 and s1a2-s2a2
        if os.path.exists(s1_allele1) and os.path.exists(s2_allele1):
            dv_s1a1_s2a1 = generate_bam_and_get_divergency(s1_allele1, s2_allele1, os.path.join(bam_output_folder, f'{sample1}_vs_{sample2}_{gene}_s1a1_s2a1.bam'))
        else:
            dv_s1a1_s2a1 = 0
        
        if os.path.exists(s1_allele2) and os.path.exists(s2_allele2):
            dv_s1a2_s2a2 = generate_bam_and_get_divergency(s1_allele2, s2_allele2, os.path.join(bam_output_folder, f'{sample1}_vs_{sample2}_{gene}_s1a2_s2a2.bam'))
        else:
            dv_s1a2_s2a2 = 0
        if dv_s1a1_s2a1 == None:
            dv_s1a1_s2a1 = 0
        if dv_s1a2_s2a2 == None:
            dv_s1a2_s2a2 = 0
        dv_s1_s2.append((dv_s1a1_s2a1+dv_s1a2_s2a2)/2)
        # for s1a1-s2a2 and s1a2-s2a1
        if os.path.exists(s1_allele1) and os.path.exists(s2_allele2):
            dv_s1a1_s2a2 = generate_bam_and_get_divergency(s1_allele1, s2_allele2, os.path.join(bam_output_folder, f'{sample1}_vs_{sample2}_{gene}_s1a1_s2a2.bam'))
        else:
            dv_s1a1_s2a2 = 0
        if os.path.exists(s1_allele2) and os.path.exists(s2_allele1):
            dv_s1a2_s2a1 = generate_bam_and_get_divergency(s1_allele2, s2_allele1, os.path.join(bam_output_folder, f'{sample1}_vs_{sample2}_{gene}_s1a2_s2a1.bam'))
        else:
            dv_s1a2_s2a1 = 0
        if dv_s1a1_s2a2 == None:
            dv_s1a1_s2a2 = 0
        if dv_s1a2_s2a1 == None:
            dv_s1a2_s2a1 = 0
        dv_s1_s2.append((dv_s1a1_s2a2+dv_s1a2_s2a1)/2)
        divergency_vector.append(min(dv_s1_s2))
    
    return divergency_vector

# 计算所有样本之间的divergency矩阵并保存为 CSV 文件，每个条目是一个divergency向量
def compute_distance_matrix(sample_names, output_csv, bam_output_folder):
    num_samples = len(sample_names)
    num_genes = len(GENES)
    
    # 初始化一个空矩阵，存储每对样本间的divergency向量
    distance_matrix = np.zeros((num_samples, num_samples, num_genes))
    
    # 确保 BAM 输出目录存在
    os.makedirs(bam_output_folder, exist_ok=True)
    
    # 遍历样本对并计算两两样本之间的divergency向量
    total_pairs = num_samples * (num_samples - 1) // 2
    for i, item in enumerate(range(num_samples)):
        for j in range(i + 1, num_samples):
            # 显示进度条
            tqdm.write(f"Processing {sample_names[i]} vs {sample_names[j]}... ({i * num_samples + j + 1}/{total_pairs})")
            divergency_vector = compute_sample_pair_divergency_vector(sample_names[i], sample_names[j], bam_output_folder)
            distance_matrix[i, j, :] = divergency_vector
            distance_matrix[j, i, :] = divergency_vector  # 距离矩阵是对称的
    
    # 展平矩阵中的每对样本的divergency向量，并保存到CSV
    with open(output_csv, 'w') as f:
        # 写入表头
        header = "Sample1,Sample2," + ",".join(GENES) + "\n"
        f.write(header)
        
        # 写入每对样本的divergency向量，避免重复输出对称样本对
        for i in range(num_samples):
            for j in range(i + 1, num_samples):
                sample1 = sample_names[i]
                sample2 = sample_names[j]
                divergency_vector = distance_matrix[i, j, :]
                line = f"{sample1},{sample2}," + ",".join(map(str, divergency_vector)) + "\n"
                f.write(line)

    print(f"Distance matrix saved to {output_csv}")

# 解析命令行参数
def parse_arguments():
    parser = argparse.ArgumentParser(description="Calculate distance matrix for HLA gene alleles between samples.")
    
    parser.add_argument('-i', '--input', required=True, help="Input file containing sample names (one per line).")
    parser.add_argument('-p', '--path', required=True, help="Base path where sample folders are stored.")
    parser.add_argument('-o', '--output', required=True, help="Output CSV file for the distance matrix.")
    parser.add_argument('-b', '--bam_output', required=True, help="Directory to save generated BAM files.")
    parser.add_argument('-t', '--threads', type=int, default=1, help="Number of threads to use for minimap2.")
    
    return parser.parse_args()

# 主函数
if __name__ == '__main__':
    # 解析命令行参数
    args = parse_arguments()
    
    # 从样本名文件读取样本名列表
    with open(args.input, 'r') as f:
        sample_names = [line.strip() for line in f if line.strip()]
    
    # 计算两两样本之间的divergency矩阵
    compute_distance_matrix(sample_names, args.output, args.bam_output)