import os
import pandas as pd
from collections import defaultdict
from tqdm import tqdm
import subprocess


gene_list = [
    'HLA-A', 'HLA-B', 'HLA-C', 'HLA-DMA', 'HLA-DMB', 'HLA-DOA', 'HLA-DOB',
    'HLA-DPA1', 'HLA-DPB1', 'HLA-DPB2', 'HLA-DQA1', 'HLA-DQB1', 'HLA-DRA',
    'HLA-DRB1', 'HLA-DRB3', 'HLA-DRB4', 'HLA-DRB5', 'HLA-E', 'HLA-F', 'HLA-G',
    'HLA-H', 'HLA-J', 'HLA-K', 'HLA-L', 'HLA-P', 'HLA-V', 'HLA-DQA2',
    'HLA-DPA2', 'HLA-N', 'HLA-S', 'HLA-T', 'HLA-U', 'HLA-W', 'MICA', 'MICB',
    'TAP1', 'TAP2', 'HFE'
]

def calculate_average_depth(bam_file):
    """
    Calculate the average depth for a BAM file using samtools and awk.

    Parameters:
    - bam_file (str): Path to the BAM file.

    Returns:
    - float: The average depth of the BAM file.
    """
    # Construct the command
    command = f"samtools depth {bam_file} | awk '{{sum+=$3}} END {{print sum/NR}}'"
    print(command)
    # Run the command
    result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, text=True)
    
    # Parse and return the output
    try:
        average_depth = float(result.stdout.strip())
    except ValueError:
        average_depth = None  # Handle cases where there's an issue with parsing the output
    
    return average_depth

def calculate_average_depth(bam_file):
    """
    Calculate the average depth for a BAM file using samtools and awk.
    Handles cases where there is no coverage.
    
    Parameters:
    - bam_file (str): Path to the BAM file.
    
    Returns:
    - float or str: The average depth of the BAM file, or a message if no coverage is found.
    """
    # Construct the command
    # '{sum+=$3} END { if (NR>0) print sum/NR; else print 0; }'
    command = f"samtools depth {bam_file} | awk '{{sum+=$3}} END {{ if (NR>0) print sum/NR; else print 0; }}'"
    
    # Run the command
    result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, text=True)
    
    # Parse the output
    output = result.stdout.strip()
    
    # Check if output is "No coverage"
    if output == "No coverage":
        return output
    
    # Convert the output to a float
    try:
        average_depth = float(output)
    except ValueError:
        average_depth = None  # Handle any unexpected output
    
    return average_depth

def count_reads_in_fastq(file_path):
    read_count = 0
    with open(file_path, 'r') as f:
        for i, line in enumerate(f):
            if i % 4 == 0:
                read_count += 1
    return read_count

def process_samples(samples_file, fastq_dirs):
    with open(samples_file, 'r') as f:
        samples = f.read().splitlines()

    dp_counts = defaultdict(lambda: defaultdict(int))

    for sample in tqdm(samples):
        for fastq_dir in fastq_dirs:
            sample_bam_dir = os.path.join(fastq_dir, sample, sample, "Genes_step2")
            if not os.path.exists(sample_bam_dir):
                continue
            for gene in gene_list:
                bam_file = os.path.join(sample_bam_dir, f'{gene}.bam')
                if os.path.exists(bam_file):
                    average_dp=calculate_average_depth(bam_file)
                    dp_counts[sample][gene] = average_dp
                else:
                    bam1_file = os.path.join(sample_bam_dir, f'{gene}.0.bam')
                    bam2_file = os.path.join(sample_bam_dir, f'{gene}.1.bam')
                    if os.path.exists(bam1_file) and os.path.exists(bam2_file):
                        average_dp1=calculate_average_depth(bam1_file)
                        average_dp2=calculate_average_depth(bam2_file)
                        dp_counts[sample][gene] = (average_dp1+average_dp2)
                    else:
                        print(f"Error: {bam_file} not found.")
                        dp_counts[sample][gene] = 0

    return dp_counts

def create_count_table(read_counts):
    df = pd.DataFrame(read_counts).fillna(0).astype(int).T
    return df

def main():
    samples_file = 'hprc1_2.samples'
    # hifi
    HPRC_hifi = "/gpfs1/scratch/ResearchGroups/cs_shuaicli/wxd/HPRC_part1/data/speclong_out_rerun/hifi"
    HPRC_ont = "/gpfs1/scratch/ResearchGroups/cs_shuaicli/wxd/HPRC_part1/data/speclong_out_rerun/ont"


    fastq_dirs = [HPRC_ont]

    dp_counts = process_samples(samples_file, fastq_dirs)
    count_table = create_count_table(dp_counts)

    print(count_table)

    # count_table.to_csv('read_counts.csv')
    count_table.to_csv('read_depth_ont.csv')

if __name__ == '__main__':
    main()