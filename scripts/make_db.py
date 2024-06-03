import os
import argparse
import subprocess
from Bio import SeqIO
from determine_gene import get_focus_gene_from_class

def download_fasta_file(url, output_path):
    """
    Download a FASTA file from a given URL and save it to a specified path using wget.

    :param url: URL of the FASTA file to download
    :param output_path: Path to save the downloaded FASTA file
    """
    result = subprocess.run(['wget', '-O', output_path, url], check=True)
    if result.returncode == 0:
        print(f"FASTA file downloaded and saved to {output_path}")
    else:
        raise Exception(f"Failed to download file from {url}")

def create_HLA_directories_and_save_sequences(fasta_path, output_base_dir, gene_list, interval_dict):
    """
    Parse a FASTA file, create directories based on gene names, 
    and save corresponding sequences to these directories.

    :param fasta_path: Path to the input FASTA file
    :param output_base_dir: Base directory to save gene-specific directories and sequences
    """
    sequences = SeqIO.parse(fasta_path, "fasta")
    gene_sequences = {}

    for seq_record in sequences:
        # Extract the gene name from the description, assuming the format is ">HLA:HLA00001 A*01:01:01:01 3503 bp"
        description_parts = seq_record.description.split()
        gene_name = description_parts[1].split('*')[0]
        gene_name = gene_name if gene_name in gene_list else f"HLA-{gene_name}"
        sequence_name = description_parts[1]

        # Update the sequence ID and name
        seq_record.id = sequence_name if gene_name in gene_list else f"HLA-{sequence_name}"
        seq_record.name = sequence_name if gene_name in gene_list else f"HLA-{sequence_name}"
        seq_record.description = ""
        
        if gene_name not in gene_sequences:
            gene_sequences[gene_name] = []
        gene_sequences[gene_name].append(seq_record)

    # Create directories and save gene-specific sequences
    all_sequences = []
    for gene_name, seq_records in gene_sequences.items():
        gene_dir = os.path.join(output_base_dir, gene_name)
        os.makedirs(gene_dir, exist_ok=True)
        gene_fasta_filename = os.path.join(gene_dir, f"{gene_name}.fasta")
        with open(gene_fasta_filename, "w") as gene_fasta_file:
            SeqIO.write(seq_records, gene_fasta_file, "fasta")
        all_sequences.extend(seq_records)
        
        # Build the index for the gene-specific FASTA file
        cmd = f"""
        samtools faidx "{gene_fasta_filename}"
        bwa index "{gene_fasta_filename}"
        {script_dir}/../bin/faToTwoBit "{gene_dir}"/{gene_name}.fasta "{gene_dir}"/{gene_name}.2bit
        makeblastdb -in "{gene_dir}"/{gene_name}.fasta -dbtype nucl -parse_seqids -out "{gene_dir}"/{gene_name}
        """
        print(cmd, flush=True)
        os.system(cmd)

    # Merge all sequences into one file and build the index
    merged_fasta_filename = os.path.join(output_base_dir, "HLA.full.fasta")
    with open(merged_fasta_filename, "w") as merged_fasta_file:
        SeqIO.write(all_sequences, merged_fasta_file, "fasta")
    
    # Build the index for the merged FASTA file
    cmd = f"""
    samtools faidx "{merged_fasta_filename}"
    bwa index "{merged_fasta_filename}"
    {script_dir}/../bin/faToTwoBit "{merged_fasta_filename}" "{os.path.join(output_base_dir, 'HLA.full.2bit')}"
    makeblastdb -in "{merged_fasta_filename}" -dbtype nucl -parse_seqids -out "{os.path.join(output_base_dir, 'HLA.full')}"
    """
    print(cmd, flush=True)
    os.system(cmd)

    print(f"All sequences have been categorized by gene and saved in {output_base_dir}")
    print(f"Merged FASTA file created and indexed at {merged_fasta_filename}")

def make_HLA_db():
    # URL to download the FASTA file
    HLA_fasta_url = "https://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/hla_gen.fasta"

    # Path to save the downloaded FASTA file within the output directory
    HLA_dir = os.path.join(args.outdir, "whole")
    print(HLA_dir)
    if not os.path.exists(HLA_dir):
        os.makedirs(HLA_dir)

    gene_list, interval_dict = get_focus_gene_from_class("HLA")
    if not args.HLA_fa:
        local_fasta_filename = os.path.join(HLA_dir, "hla_gen.fasta")
        download_fasta_file(HLA_fasta_url, local_fasta_filename)
    else:
        local_fasta_filename = args.HLA_fa

    # Parse the FASTA file and save sequences by gene
    create_HLA_directories_and_save_sequences(local_fasta_filename, HLA_dir, gene_list, interval_dict)

def main():
    make_HLA_db()

if __name__ == "__main__":
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Download and process HLA gene sequences.")
    required = parser.add_argument_group("Required arguments")
    optional = parser.add_argument_group("Optional arguments")
    
    required.add_argument("-o","--outdir", help="Directory to save the gene-specific sequences")
    # add default hla_gene.fa
    optional.add_argument("--HLA_fa", help="hla_gene")
    args = parser.parse_args()
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Ensure the output directory exists and is not a file
    if os.path.isfile(args.outdir):
        print(f"The specified output directory '{args.outdir}' is a file. Deleting the file.")
        os.remove(args.outdir)

    if args.outdir and not os.path.exists(args.outdir):
        print("No DB dir, creating one")
        os.makedirs(args.outdir)

    main()