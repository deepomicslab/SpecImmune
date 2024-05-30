import os
import pysam

def process_fasta_files(input_folder, output_folder):
    # Create the output folder if it doesn't exist
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Get the list of input files in the input folder
    input_files = [f for f in os.listdir(input_folder) if os.path.isfile(os.path.join(input_folder, f))]

    # Process each input file
    for input_file in input_files:

        # Check if the input file has a .fasta suffix
        if not input_file.endswith('.fasta'):
            continue

        # Construct the input and output file paths
        input_file_path = os.path.join(input_folder, input_file)
        output_file_path = os.path.join(output_folder, input_file)

        # Process the input file and save the modified sequences to the output file
        with pysam.FastaFile(input_file_path) as input_fasta, open(output_file_path, 'w') as output_fasta:
            contig_ids = input_fasta.references

            for contig_id in contig_ids:
                sequence = input_fasta.fetch(contig_id)
                modified_sequence = sequence[300:-300]

                output_fasta.write(f">{contig_id}\n{modified_sequence}\n")


db = "/mnt/d/HLAPro_backup/Nanopore_optimize/SpecComplex/db/"

# Usage example
input_folder = db + 'whole/'  # Replace with the actual input folder path
output_folder = db + 'clean_whole/'  # Replace with the desired output folder name


process_fasta_files(input_folder, output_folder)