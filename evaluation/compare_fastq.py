import gzip
import pysam

def count_shared_reads(file1, file2):
    # Read the gzipped FASTQ files and store the read names
    with gzip.open(file1, 'rt') as gzipped_file1, gzip.open(file2, 'rt') as gzipped_file2:
        file1_path = gzipped_file1.name
        file2_path = gzipped_file2.name

        fastq1 = pysam.FastxFile(file1_path)
        fastq2 = pysam.FastxFile(file2_path)

        read_names1 = set(record.name for record in fastq1)
        read_names2 = set(record.name for record in fastq2)

    # Calculate the number of shared reads
    shared_reads = read_names1.intersection(read_names2)

    # Calculate the file-specific reads
    file1_specific_reads = read_names1.difference(read_names2)
    file2_specific_reads = read_names2.difference(read_names1)

    return shared_reads, file1_specific_reads, file2_specific_reads

# Usage example
gene="DQB1"
gene_list = ['A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRB1']
for gene in gene_list:
    print (gene)
    file1 = f'/mnt/d/HLAPro_backup/Nanopore_optimize/output4/fredhutch-hla-FH5_0.1/HLA-{gene}.long_read.fq.gz'  # Replace with the actual gzipped file path of the first FASTQ file
    file2 = f'/mnt/d/HLAPro_backup/Nanopore_optimize/output0/fredhutch-hla-FH5_0.1/{gene}.long_read.fq.gz'  # Replace with the actual gzipped file path of the second FASTQ file

    shared_reads, file1_specific_reads, file2_specific_reads = count_shared_reads(file1, file2)

    print("Number of shared reads:", len(shared_reads))
    print("Number of file1-specific reads:", len(file1_specific_reads))
    print("Number of file2-specific reads:", len(file2_specific_reads))


    # Print the file-specific read names
    print("\nFile1-specific reads:")
    for name in file1_specific_reads:
        print(name)

    print("\nFile2-specific reads:")
    for name in file2_specific_reads:
        print(name)