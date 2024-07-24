import requests
import os
import argparse
from tqdm import tqdm
import sys

# Base URL for the FTP server
base_url = 'https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/'

# Function to download a file with progress bar
def download_file(file_url, file_path, file_size):
    response = requests.get(file_url, stream=True)
    if response.status_code == 200:
        with open(file_path, 'wb') as f, tqdm(
            desc=file_path,
            total=int(file_size),
            unit='B',
            unit_scale=True,
            unit_divisor=1024,
        ) as bar:
            for chunk in response.iter_content(1024):
                if chunk:  # filter out keep-alive new chunks
                    f.write(chunk)
                    bar.update(len(chunk))
        print(f"Downloaded {file_url} to {file_path}")
    else:
        print(f"Failed to download {file_url}. Status code: {response.status_code}")

# Function to get the sample name from the file path
def get_sample_name(file_path):
    # Extract the filename
    filename = os.path.basename(file_path)
    # Assuming the sample name is the part before the first period
    sample_name = filename.split('.')[0]
    return sample_name

def load_GeT_RM(given_sample_list):
    GeT_RM_sample_list = []
    sampl_truth = {}
    with open('cyp/CYP_GeT-RM_truth.csv') as f:
        for line in f:
            line = line.strip().split(',')
            sample = line[0]
            if len(line) < 4:
                continue
            truth = line[3]
            sampl_truth[sample] = truth
            # print (sample, truth)
            if sample in given_sample_list:
                GeT_RM_sample_list.append(sample)
    print (len(GeT_RM_sample_list))
    print (GeT_RM_sample_list)
    # file = open('cyp/ont_truth.csv', 'w')
    # for sample in GeT_RM_sample_list:
    #     file.write(f"{sample}\t{sampl_truth[sample]}\n")
    # file.close()

def load_GeT_RM3(given_sample_list):
    GeT_RM_sample_list = []
    sampl_truth = {}
    with open('cyp/GeT_RM_truth.txt') as f:
        for line in f:
            line = line.strip().split('\t')
            sample = line[0]
            print (sample, line[18])
            if len(line) < 4:
                continue
            truth = line[17]
            sampl_truth[sample] = truth
            # print (sample, truth)
            if sample in given_sample_list:
                GeT_RM_sample_list.append(sample)
    print (len(GeT_RM_sample_list))
    print (GeT_RM_sample_list)
    # file = open('cyp/ont_truth.csv', 'w')
    # for sample in GeT_RM_sample_list:
    #     file.write(f"{sample}\t{sampl_truth[sample]}\n")
    # file.close()


def load_GeT_RM2(given_sample_list):
    GeT_RM_sample_list = []
    with open('cyp/cyp_benchmark.csv') as f:
        for line in f:
            line = line.strip().split()
            sample = line[1]
            # truth = line[1]
            # print (sample, truth)
            if sample in given_sample_list:
                GeT_RM_sample_list.append(sample)
    print (len(GeT_RM_sample_list))
    print (GeT_RM_sample_list)

# Main function to process the file list and download files
def main(file_list_path):
    # Read the file list
    ont_sample_list = []
    with open(file_list_path, 'r') as file:
        for line in file:
            # Each line in the text file is expected to have the format: path size checksum
            parts = line.strip().split('\t')
            if len(parts) != 3:
                print(f"Skipping invalid line: {line}")
                continue

            file_path, file_size, checksum = parts

            # Only download files with the 'hg38' tag
            if 'hg38' in file_path:
                file_url = base_url + file_path
                sample_name = get_sample_name(file_path)
                ont_sample_list.append(sample_name)
                # print(f"Downloading {file_path} for sample {sample_name}")
                for gene_class in ['HLA', 'KIR', 'CYP', 'IG_TR']:
                    cmd = f"""
                        bash {sys.path[0]}/../scripts/ExtractReads.sh -s {sample_name} -i /mnt/e/1000G_ONT/1000G_ont/downloads/{sample_name}/{sample_name}.hg38.cram -g {gene_class} -o /mnt/d/HLAPro_backup/Nanopore_optimize/data/1000G_ont/downloads/{sample_name} -r /mnt/d/HLAPro_backup/Nanopore_optimize/data/hg38/1KG_ONT_VIENNA_hg38.fa
                    """
                    if not os.path.exists(f"/mnt/e/1000G_ONT/1000G_ont/downloads/{sample_name}/{sample_name}.{gene_class}.fastq.gz"):
                        # print (cmd)
                        os.system(cmd)
    # load_GeT_RM(ont_sample_list)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Download files from a given list.')
    parser.add_argument('file_list', type=str, help='Path to the file containing the list of files to download.')
    args = parser.parse_args()

    main(args.file_list)
