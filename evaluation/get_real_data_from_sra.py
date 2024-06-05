import os
import csv

# Function to extract ERR number from URL
def extract_err_number(url):
    print(url)
    parts = url.split('/')
    print(parts)
    return parts[6]

# Function to extract the original BAM file name from the URL
def extract_bam_filename(url):
    parts = url.split('/')
    return parts[-1]

# Function to extract the sample name from the BAM file name
def extract_sample_name(filename):
    parts = filename.replace("-", "_").split("_")
    return parts[0]

# Function to download file using prefetch, convert to FASTQ format using fasterq-dump
def download_and_convert(url, threads=4):
    err_number = extract_err_number(url)
    bam_filename = extract_bam_filename(url)
    sample_name = extract_sample_name(bam_filename)
    print(err_number, bam_filename, sample_name)

    # Create a directory for the sample
    if not os.path.exists(sample_name):
        os.makedirs(sample_name)
    
    # Download the file using prefetch
    prefetch_command = f"prefetch {err_number} --max-size 800G -p -O {sample_name}"
    os.system(prefetch_command)

    # Convert the downloaded SRA file to FASTQ format using fasterq-dump
    fasterq_dump_command = f"fasterq-dump --threads {threads} {sample_name}/{err_number}/{err_number}.sra -O {sample_name}/{err_number}"
    os.system(fasterq_dump_command)

    # Remove the SRA file to save space
    # os.remove(f"{sample_name}/{err_number}/{err_number}.sra")

# Input file path
input_file = 'igsr_HGSVC2_PacBio_CLR.txt'

# Read the CSV file
with open(input_file, 'r') as file:
    reader = csv.reader(file, delimiter='\t')
    headers = next(reader)  # Skip the header row

    for row in reader:
        url = row[0]
        download_and_convert(url)