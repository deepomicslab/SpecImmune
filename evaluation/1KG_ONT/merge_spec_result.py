import os
import pandas as pd
from tqdm import tqdm

def process_samples(samples_file, fastq_dirs):
    """
    Processes a list of samples, reads their corresponding result files, and 
    merges them into a single DataFrame, with an additional 'Sample' column.

    Args:
        samples_file (str): Path to the file containing the list of sample names.
        fastq_dirs (list): List of directories where the sample results are stored.

    Returns:
        pd.DataFrame: DataFrame containing all merged sample results with a 'Sample' column.
    """
    # Read all sample names from the samples file
    with open(samples_file, 'r') as f:
        samples = f.read().splitlines()

    all_samples_df = pd.DataFrame()

    # Iterate over each sample and directory
    for sample in tqdm(samples, desc="Processing samples"):
        for fastq_dir in fastq_dirs:
            sample_res_dir = os.path.join(fastq_dir, sample, sample)
            sample_res_file = os.path.join(sample_res_dir, f"{sample}.HLA.final.type.result.formatted.txt")

            if os.path.exists(sample_res_file):
                try:
                    # Read the sample result file, skipping the first line (header comment)
                    df = pd.read_csv(sample_res_file, sep='\t', skiprows=1)
                    df['Sample'] = sample  # Add a column for the sample name
                    all_samples_df = pd.concat([all_samples_df, df], ignore_index=True)
                except Exception as e:
                    print(f"Error reading {sample_res_file}: {e}")
            else:
                print(f"File {sample_res_file} does not exist.")

    return all_samples_df

def main():
    """
    Main function to execute the processing of samples and generation of count tables.
    """
    # Define the samples file and result directories
    samples_file = '3parts.merge.samples.unique'
    KGP_ONT_part1 = "/gpfs1/scratch/ResearchGroups/cs_shuaicli/wxd/1KGP_ONT/speclong_out_rerun/"
    # Additional directories can be added here
    res_dirs = [KGP_ONT_part1]

    # Process samples and create the merged DataFrame
    all_samples_df = process_samples(samples_file, res_dirs)



    # Save the merged DataFrame and count table to CSV files
    all_samples_df.to_csv('merged_samples.csv', index=False)

if __name__ == '__main__':
    main()