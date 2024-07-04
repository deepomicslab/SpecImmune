import requests
import os
import argparse
from tqdm import tqdm

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

# Main function to process the file list and download files
def main(file_list_path):
    # Read the file list
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

                # Directory to save the downloaded files for this sample
                sample_dir = os.path.join('downloads', sample_name)
                os.makedirs(sample_dir, exist_ok=True)

                local_file_path = os.path.join(sample_dir, os.path.basename(file_path))

                # Download the file
                download_file(file_url, local_file_path, file_size)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Download files from a given list.')
    parser.add_argument('file_list', type=str, help='Path to the file containing the list of files to download.')
    args = parser.parse_args()

    main(args.file_list)
