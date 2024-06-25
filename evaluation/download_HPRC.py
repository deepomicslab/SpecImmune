import csv
import os
import subprocess
from tqdm import tqdm

def get_file_size(s3_url):
    """Get the size of a file from an S3 URL using aws s3 ls command with --no-sign-request."""
    try:
        command = ['aws', 's3', 'ls', s3_url, '--no-sign-request']
        result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        output = result.stdout.strip()
        if result.returncode != 0:
            print(f"Failed to get file size for {s3_url}: {result.stderr}")
            return None
        if output:
            size = int(output.split()[2])
            return size
        return None
    except subprocess.CalledProcessError as e:
        print(f"Failed to get file size for {s3_url}: {e}")
        return None

def download_from_s3(s3_url, local_path):
    """Download a file from an S3 URL to a local path using aws s3 cp command with --no-sign-request."""
    try:
        command = ['aws', 's3', 'cp', s3_url, local_path, '--no-sign-request']
        subprocess.run(command, check=True)
        print(f"Downloaded {s3_url} to {local_path}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"Failed to download {s3_url}: {e}")
        return False

def download_files(urls, sample_id, data_type, max_size_gb=300):
    """Download files from S3 based on size constraints."""
    downloaded_size = 0
    max_size_bytes = max_size_gb * 1024 ** 3
    failed_downloads = []

    for i, url in enumerate(tqdm(urls, desc=f"Downloading {data_type} files for {sample_id}")):
        if not url:
            continue

        size = get_file_size(url)
        if size is None:
            failed_downloads.append(url)
            continue

        if downloaded_size + size > max_size_bytes:
            break

        local_path = os.path.join('data', sample_id, data_type, os.path.basename(url))
        os.makedirs(os.path.dirname(local_path), exist_ok=True)
        if download_from_s3(url, local_path):
            downloaded_size += size
            print(f"Downloaded {url} ({size} bytes), total downloaded size: {downloaded_size} bytes")
        else:
            failed_downloads.append(url)

    return failed_downloads

def process_tsv(tsv_file):
    """Process the TSV file and download HiFi and Nanopore data."""
    all_failed_downloads = []

    with open(tsv_file, newline='', encoding='utf-8-sig') as file:  # Use utf-8-sig to handle BOM
        reader = csv.DictReader(file, delimiter=',')
        headers = reader.fieldnames
        print(f"Headers found in TSV: {headers}")

        if 'sample_id' not in headers:
            raise KeyError("The TSV file does not contain a 'sample_id' column.")

        for row in reader:
            sample_id = row['sample_id']
            hifi_urls = row['hifi'].split(';') if row['hifi'] else []
            nanopore_urls = row['nanopore'].split(';') if row['nanopore'] else []

            # Download HiFi data
            failed_downloads = download_files(hifi_urls, sample_id, 'hifi')
            all_failed_downloads.extend(failed_downloads)

            # Download Nanopore data
            failed_downloads = download_files(nanopore_urls, sample_id, 'nanopore')
            all_failed_downloads.extend(failed_downloads)

    # Print failed downloads
    if all_failed_downloads:
        print("\nFailed downloads:")
        for url in all_failed_downloads:
            print(url)
    else:
        print("\nAll downloads completed successfully.")

if __name__ == "__main__":
    tsv_file = 'forth_part.csv'  # Replace with your TSV file path
    process_tsv(tsv_file)
