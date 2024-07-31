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

def load_GeT_RM4():
    GeT_RM_truth = {}
    ## read the excel into dataframe
    import pandas as pd
    file = "cyp/GeT_RM_truth.csv"
    out_file = "cyp/GeT_RM_truth.tab"
    f = open(out_file, 'w')
    df = pd.read_csv(file)
    # print (df.columns)
    ## enumerate each row
    for index, row in df.iterrows():
        sample = row['Coriell #  https://www.coriell.org/  ']
        truth = row['CYP2D6']
        ## check if truth is nan
        if pd.isna(truth):
            continue
        if truth == 'no consensus' or truth == 'duplication':
            continue
        print (sample, truth, file = f)
        GeT_RM_truth[sample] = truth
    f.close()
    return GeT_RM_truth

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
                if sample_name not in ont_sample_list:
                    ont_sample_list.append(sample_name)
                # print(f"Downloading {file_path} for sample {sample_name}")
                # for gene_class in ['HLA', 'KIR', 'CYP', 'IG_TR']:
                #     cmd = f"""
                #         bash {sys.path[0]}/../scripts/ExtractReads.sh -s {sample_name} -i /mnt/e/1000G_ONT/1000G_ont/downloads/{sample_name}/{sample_name}.hg38.cram -g {gene_class} -o /mnt/d/HLAPro_backup/Nanopore_optimize/data/1000G_ont/downloads/{sample_name} -r /mnt/d/HLAPro_backup/Nanopore_optimize/data/hg38/1KG_ONT_VIENNA_hg38.fa
                #     """
                #     if not os.path.exists(f"/mnt/e/1000G_ONT/1000G_ont/downloads/{sample_name}/{sample_name}.{gene_class}.fastq.gz"):
                #         # print (cmd)
                #         os.system(cmd)
    # load_GeT_RM(ont_sample_list)
    return ont_sample_list

def check_share_sample(GeT_RM_truth, ont_sample_list, file):
    file = open(file, 'w')
    # for sample in GeT_RM_sample_list:
    shared_num = 0
    for sample in ont_sample_list:
        if sample in GeT_RM_truth:
            print (sample, GeT_RM_truth[sample])
            file.write(f"{sample}\t{GeT_RM_truth[sample]}\n")
            shared_num += 1
    file.close()
    # print (ont_sample_list)
    print (shared_num)

def load_hgscv_sample(file):
    sample_list = []
    ## load file by pd
    import pandas as pd
    df = pd.read_csv(file, sep = '\t')
    ## for  each row
    for index, row in df.iterrows():
        sample = row['Sample']
        if sample not in sample_list:
            sample_list.append(sample)
        # print (sample)
    print ("sample num", len(sample_list))
    return sample_list

def load_hprc_sample(file):
    sample_list = []
    ## load file by pd
    import pandas as pd
    df = pd.read_csv(file, sep = ',')
    ## for  each row
    for index, row in df.iterrows():
        sample = row['sample_id']
        if sample not in sample_list:
            sample_list.append(sample)
        # print (sample)
    print ("sample num", len(sample_list))
    print (sample_list)
    return sample_list

def load_Stargazer_truth():
    Stargazer_truth = {}
    ## load file by pd
    import pandas as pd
    df = pd.read_csv('cyp/stargazer.csv', sep = ',')
    ## for  each row
    for index, row in df.iterrows():
        # print (index, row)
        sample = row['Sample']
        truth = row['Orthogonal methods']
        Stargazer_truth[sample] = truth
        # print (sample)
    # print ("sample num", len(Stargazer_truth))

    return Stargazer_truth


if __name__ == '__main__':
    # parser = argparse.ArgumentParser(description='Download files from a given list.')
    # parser.add_argument('file_list', type=str, help='Path to the file containing the list of files to download.')
    # args = parser.parse_args()
    
    GeT_RM_truth1 = load_GeT_RM4()
    GeT_RM_truth2 = load_Stargazer_truth()

    GeT_RM_truth = {**GeT_RM_truth1, **GeT_RM_truth2}
    ont_file_list = "lkg_ont_vienna_merge.files.hg38.list"
    ont_sample_list = main(ont_file_list)
    truth_file = 'cyp/ont_truth_merge.csv'
    check_share_sample(GeT_RM_truth, ont_sample_list, truth_file)

    # hgscv_clr = "igsr_HGSVC2_PacBio_CLR.tsv"
    # sample_list = load_hgscv_sample(hgscv_clr)
    # truth_file = 'cyp/HGSVC2_PacBio_CLR_truth.csv'
    # check_share_sample(GeT_RM_truth, sample_list, truth_file)

    # hgscv_clr = "igsr_HGSVC2_PacBio_HIFI.tsv"
    # sample_list = load_hgscv_sample(hgscv_clr)
    # truth_file = 'cyp/HGSVC2_PacBio_HIFI_truth.csv'
    # check_share_sample(GeT_RM_truth, sample_list, truth_file)

    # hgscv_clr = "hprc_year1_level1_aws_locs.csv"
    # sample_list = load_hprc_sample(hgscv_clr)
    # truth_file = 'cyp/hprc_year1_level1_truth.csv'
    # check_share_sample(GeT_RM_truth, sample_list, truth_file)
