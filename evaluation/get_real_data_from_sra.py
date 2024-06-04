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

# Function to download file using prefetch, convert to SAM, and then to BAM format
def download_and_convert(url):
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

    # Convert the downloaded SRA file to SAM format using sam-dump
    sam_dump_command = f"sam-dump {sample_name}/{err_number}/{err_number}.sra > {sample_name}/{err_number}/{err_number}.sam"
    os.system(sam_dump_command)
    
    # Convert SAM to BAM format using samtools and move to the sample directory
    bam_command = f"samtools view -Sb {sample_name}/{err_number}/{err_number}.sam > {sample_name}/{err_number}/{bam_filename}"
    os.system(bam_command)

    # Optionally, remove the intermediate SAM file to save space
    os.remove(f"{sample_name}/{err_number}/{err_number}.sam")

# Input file path
input_file = 'igsr_HGSVC2_PacBio_CLR.txt'

# Read the CSV file
with open(input_file, 'r') as file:
    reader = csv.reader(file, delimiter='\t')
    headers = next(reader)  # Skip the header row

    for row in reader:
        url = row[0]
        download_and_convert(url)