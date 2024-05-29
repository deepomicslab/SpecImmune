import urllib.request
import hashlib
import os

def download_and_check_md5(url, md5, data_dir):
    # Extract the filename from the URL
    filename = url.split('/')[-1]
    sample = filename[:-4]
    # Download the file
    # urllib.request.urlretrieve(url, os.path.join(data_dir, filename))
    os.system(f"wget {url} -O {os.path.join(data_dir, filename)}")
    

    # Calculate MD5 checksum of the downloaded file
    with open(os.path.join(data_dir, filename), 'rb') as file:
        file_content = file.read()
        calculated_md5 = hashlib.md5(file_content).hexdigest()

    # Compare the calculated MD5 with the provided MD5
    if calculated_md5 == md5:
        print(f"MD5 checksum matches for {filename}")
        bam = os.path.join(data_dir, filename)
        cmd = f"""
        samtools index {bam}
        pbmm align {hg38} {bam} {data_dir}/{sample}.hg38.bam --sort -j 16 -J 8
        bash ../scripts/run.extract.reads.sh {sample} {data_dir}/{sample}.hg38.bam HLA {data_dir}
        bash ../scripts/run.extract.reads.sh {sample} {data_dir}/{sample}.hg38.bam KIR {data_dir}
        bash ../scripts/run.extract.reads.sh {sample} {data_dir}/{sample}.hg38.bam CYP {data_dir}
        bash ../scripts/run.extract.reads.sh {sample} {data_dir}/{sample}.hg38.bam IG_TR {data_dir}
        """
    else:
        print(f"MD5 checksum does not match for {filename}")

        # Remove the downloaded file
        # Comment out the following line if you want to keep the downloaded file
        # after checking the MD5 checksum
        os.remove(os.path.join(data_dir, filename))


# Read the data list file
# data_dir = '/mnt/d/HLAPro_backup/Nanopore_optimize/data/pacbio_clr/'  # Specify the data directory here
data_dir = "/RGC_1/shuai/hla_pacbio_clr/"
data_list_file = 'igsr_HGSVC2_PacBio_CLR.tsv'  # Specify the data list file here
hg38="/home/wangmengyao/GATK/database/Homo_sapiens_assembly38.fasta"

with open(data_list_file, 'r') as file:
    file.readline()
    for line in file:
        line = line.strip()
        if line:
            # Split the line by tab or space
            data = line.split('\t') if '\t' in line else line.split()

            # Extract the URL and MD5 checksum
            url = data[0]
            md5 = data[1]

            # Download and check MD5 checksum
            download_and_check_md5(url, md5, data_dir)
