import os
import argparse
import subprocess
from Bio import SeqIO

def download_file(url, output_path):
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

def create_HLA_directories_and_save_sequences(fasta_path, output_base_dir, gene_list):
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
        # gene_name = gene_name if gene_name in gene_list else f"HLA-{gene_name}"
        sequence_name = description_parts[1]

        # Update the sequence ID and name
        seq_record.id = sequence_name if gene_name in gene_list else f"HLA-{sequence_name}"
        seq_record.name = sequence_name if gene_name in gene_list else f"HLA-{sequence_name}"
        seq_record.description = ""
        gene_name = gene_name if gene_name in gene_list else f"HLA-{gene_name}"

        print(seq_record.id, seq_record.name, seq_record.description, flush=True)
        
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
    makeblastdb -in "{merged_fasta_filename}" -dbtype nucl -parse_seqids -out "{os.path.join(output_base_dir, 'HLA.full')}"
    """
    print(cmd, flush=True)
    os.system(cmd)

    print(f"All sequences have been categorized by gene and saved in {output_base_dir}")
    print(f"Merged FASTA file created and indexed at {merged_fasta_filename}")

def create_KIR_directories_and_save_sequences(fasta_path, output_base_dir, gene_list, interval_dict):
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
        sequence_name = description_parts[1]
        

        # Update the sequence ID and name
        seq_record.id = sequence_name if gene_name in gene_list else f"{sequence_name}"
        seq_record.name = sequence_name if gene_name in gene_list else f"{sequence_name}"
        seq_record.description = ""
        # if gene_name.startswith("KIR2DL5"):
        if False:
            gene_name = gene_name
        elif gene_name in gene_list:
            gene_name = gene_name
        else:
            gene_name = f"KIR-{gene_name}"

        print(seq_record.id, seq_record.name, seq_record.description, flush=True)
        
        if gene_name not in gene_sequences:
            gene_sequences[gene_name] = []
        gene_sequences[gene_name].append(seq_record)

    # Create directories and save gene-specific sequences
    all_sequences = []
    for gene_name, seq_records in gene_sequences.items():
        gene_dir = os.path.join(output_base_dir, gene_name)
        # if gene_name.startswith("KIR2DL5"):
        #     gene_dir = os.path.join(output_base_dir, "KIR2DL5")

        os.makedirs(gene_dir, exist_ok=True)
        gene_fasta_filename = os.path.join(gene_dir, f"{gene_name}.fasta")
        with open(gene_fasta_filename, "w") as gene_fasta_file:
            SeqIO.write(seq_records, gene_fasta_file, "fasta")
        all_sequences.extend(seq_records)
        
        # Build the index for the gene-specific FASTA file
        cmd = f"""
        samtools faidx "{gene_fasta_filename}"
        bwa index "{gene_fasta_filename}"
        makeblastdb -in "{gene_dir}"/{gene_name}.fasta -dbtype nucl -parse_seqids -out "{gene_dir}"/{gene_name}
        """
        print(cmd, flush=True)
        os.system(cmd)

    # Merge all sequences into one file and build the index
    merged_fasta_filename = os.path.join(output_base_dir, "KIR.full.fasta")
    with open(merged_fasta_filename, "w") as merged_fasta_file:
        SeqIO.write(all_sequences, merged_fasta_file, "fasta")

    # Build the index for the merged FASTA file
    cmd = f"""
    samtools faidx "{merged_fasta_filename}"
    bwa index "{merged_fasta_filename}"
    makeblastdb -in "{merged_fasta_filename}" -dbtype nucl -parse_seqids -out "{os.path.join(output_base_dir, 'KIR.full')}"
    """
    print(cmd, flush=True)
    os.system(cmd)

    print(f"All sequences have been categorized by gene and saved in {output_base_dir}")
    print(f"Merged FASTA file created and indexed at {merged_fasta_filename}")

def create_CYP_directories_and_save_sequences(fasta_path, output_base_dir, gene_list):
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
        allele_name = seq_record.id
        gene_name = allele_name.split('*')[0]
        if gene_name not in gene_list:
            continue
        

        # Update the sequence ID and name
        # seq_record.id = sequence_name if gene_name in gene_list else f"{sequence_name}"
        # seq_record.name = sequence_name if gene_name in gene_list else f"{sequence_name}"
        seq_record.description = ""
        # if gene_name.startswith("KIR2DL5"):
        if False:
            gene_name = gene_name
        elif gene_name in gene_list:
            gene_name = gene_name
        else:
            gene_name = f"CYP-{gene_name}"

        print(seq_record.id, seq_record.name, seq_record.description, flush=True)
        
        if gene_name not in gene_sequences:
            gene_sequences[gene_name] = []
        gene_sequences[gene_name].append(seq_record)

    # Create directories and save gene-specific sequences
    all_sequences = []
    for gene_name, seq_records in gene_sequences.items():
        gene_dir = os.path.join(output_base_dir, gene_name)
        # if gene_name.startswith("KIR2DL5"):
        #     gene_dir = os.path.join(output_base_dir, "KIR2DL5")

        os.makedirs(gene_dir, exist_ok=True)
        gene_fasta_filename = os.path.join(gene_dir, f"{gene_name}.fasta")
        with open(gene_fasta_filename, "w") as gene_fasta_file:
            SeqIO.write(seq_records, gene_fasta_file, "fasta")
        all_sequences.extend(seq_records)
        
        # Build the index for the gene-specific FASTA file
        cmd = f"""
        samtools faidx "{gene_fasta_filename}"
        bwa index "{gene_fasta_filename}"
        makeblastdb -in "{gene_dir}"/{gene_name}.fasta -dbtype nucl -parse_seqids -out "{gene_dir}"/{gene_name}
        """
        print(cmd, flush=True)
        os.system(cmd)

    # Merge all sequences into one file and build the index
    merged_fasta_filename = os.path.join(output_base_dir, "CYP.full.fasta")
    with open(merged_fasta_filename, "w") as merged_fasta_file:
        SeqIO.write(all_sequences, merged_fasta_file, "fasta")

    # Build the index for the merged FASTA file
    cmd = f"""
    samtools faidx "{merged_fasta_filename}"
    bwa index "{merged_fasta_filename}"
    makeblastdb -in "{merged_fasta_filename}" -dbtype nucl -parse_seqids -out "{os.path.join(output_base_dir, 'CYP.full')}"
    """
    print(cmd, flush=True)
    os.system(cmd)

    print(f"All sequences have been categorized by gene and saved in {output_base_dir}")
    print(f"Merged FASTA file created and indexed at {merged_fasta_filename}")

def make_HLA_db():
    # URL to download the FASTA file
    HLA_fasta_url = "https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/hla_gen.fasta"
    release_version = "https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/release_version.txt"

    # Path to save the downloaded FASTA file within the output directory
    HLA_dir = os.path.join(args.outdir, "HLA")
    print(HLA_dir)
    if not os.path.exists(HLA_dir):
        os.makedirs(HLA_dir)

    if not args.HLA_fa:
        local_fasta_filename = os.path.join(HLA_dir, "hla_gen.fasta")
        local_release_version = os.path.join(HLA_dir, "release_version.txt")
        download_file(HLA_fasta_url, local_fasta_filename)
        download_file(release_version, local_release_version)
    else:
        local_fasta_filename = args.HLA_fa
    gene_list = []
    for record in SeqIO.parse(local_fasta_filename, "fasta"):
        ctg_name=record.id.split('*')[0]
        if ctg_name.startswith('rs'):
            continue
        gene_list.append(ctg_name)
    

    # Parse the FASTA file and save sequences by gene
    create_HLA_directories_and_save_sequences(local_fasta_filename, HLA_dir, gene_list)


def remove_duplicate_contigs(fasta_file, output_file):
    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        if str(record.seq) not in sequences:
            sequences[str(record.seq)] = record
    SeqIO.write(sequences.values(), output_file, "fasta")

def make_KIR_db():
    # URL to download the FASTA file
    KIR_fasta_url = "https://raw.githubusercontent.com/ANHIG/IPDKIR/Latest/kir_gen.fasta"
    release_version = "https://raw.githubusercontent.com/ANHIG/IPDKIR/Latest/release_version.txt"

    # Path to save the downloaded FASTA file within the output directory
    KIR_dir = os.path.join(args.outdir, "KIR")
    print(KIR_dir)
    if not os.path.exists(KIR_dir):
        os.makedirs(KIR_dir)

    if not args.KIR_fa:
        local_fasta_filename = os.path.join(KIR_dir, "kir_gen.fasta")
        local_release_version = os.path.join(KIR_dir, "release_version.txt")
        download_file(KIR_fasta_url, local_fasta_filename)
        download_file(release_version, local_release_version)
    else:
        local_fasta_filename = args.KIR_fa
    
    unique_fasta_filename = os.path.join(KIR_dir, "kir_gen_unique.fasta")
    remove_duplicate_contigs(local_fasta_filename, unique_fasta_filename)
    gene_list=[]
    for record in SeqIO.parse(unique_fasta_filename, "fasta"):
        ctg_name=record.id.split('*')[0]
        if ctg_name.startswith('rs'):
            continue
        gene_list.append(ctg_name)
    # Parse the FASTA file and save sequences by gene
    create_KIR_directories_and_save_sequences(unique_fasta_filename, KIR_dir, gene_list)

def make_CYP_db():
    # Path to save the downloaded FASTA file within the output directory
    CYP_dir = os.path.join(args.outdir, "CYP")
    print(CYP_dir)
    if not os.path.exists(CYP_dir):
        os.makedirs(CYP_dir)

    if not args.CYP_fa:
        local_fasta_filename = os.path.join(CYP_dir, "cyp_gen.fasta")
        local_release_version = os.path.join(CYP_dir, "release_version.txt")
    else:
        local_fasta_filename = args.CYP_fa
    # change gene_list to the contig names in the fasta file
    gene_list = []
    for record in SeqIO.parse(local_fasta_filename, "fasta"):
        ctg_name=record.id.split('*')[0]
        if ctg_name.startswith('rs'):
            continue
        gene_list.append(ctg_name)
    
    # Parse the FASTA file and save sequences by gene
    create_CYP_directories_and_save_sequences(local_fasta_filename, CYP_dir, gene_list)


def main():
    # make_HLA_db()
    # make_KIR_db()
    if args.HLA_fa:
        make_HLA_db()
    if args.KIR_fa:
        make_KIR_db()
    if args.CYP_fa:
        make_CYP_db()

if __name__ == "__main__":
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Download and process HLA gene sequences.")
    required = parser.add_argument_group("Required arguments")
    optional = parser.add_argument_group("Optional arguments")
    
    required.add_argument("-o","--outdir", help="Directory to save the gene-specific sequences")
    # add default hla_gene.fa
    optional.add_argument("--HLA_fa", help="hla_gene")
    optional.add_argument("--KIR_fa", help="kir_gene")
    optional.add_argument("--CYP_fa", help="cyp_gene")
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