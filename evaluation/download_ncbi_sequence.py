# from Bio import Entrez, SeqIO

# # Always provide your email when using NCBI's Entrez
# Entrez.email = "wshuai294@gmail.com"

# # The accession number of the sequence you want to download
# accession = "OP255873.1"

# # Fetch the sequence from NCBI
# with Entrez.efetch(db="nucleotide", rettype="fasta", retmode="text", id=accession) as handle:
#     # Read the fetched sequence
#     seq_record = SeqIO.read(handle, "fasta")

# # Save the sequence to a file
# output_file = f"/mnt/d/HLAPro_backup/Nanopore_optimize/data/VDJ/{accession}.fasta"
# with open(output_file, "w") as f:
#     SeqIO.write(seq_record, f, "fasta")

# print(f"Downloaded {accession} and saved as {output_file}")


from Bio import Entrez, SeqIO

# Always provide your email when using NCBI's Entrez
Entrez.email = "wshuai294@gmail.com"

# Define your search term
search_term = "555323[BioProject]"

# Use Entrez.esearch to get a list of accession numbers for your search term
with Entrez.esearch(db="nucleotide", term=search_term, retmax=100000) as search_handle:
    search_results = Entrez.read(search_handle)
    id_list = search_results["IdList"]

# Directory to save the downloaded sequences
output_dir = "/mnt/d/HLAPro_backup/Nanopore_optimize/data/VDJ/ncbi_assembly/"

# Iterate through each ID and download the corresponding sequence
for accession in id_list:
    with Entrez.efetch(db="nucleotide", rettype="fasta", retmode="text", id=accession) as handle:
        seq_record = SeqIO.read(handle, "fasta")
        output_file = f"{output_dir}/{accession}.fasta"
        with open(output_file, "w") as f:
            SeqIO.write(seq_record, f, "fasta")
        print(f"Downloaded {accession} and saved as {output_file}")