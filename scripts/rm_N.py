import argparse
from Bio import SeqIO

def remove_N_characters(input_file, output_file):
    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        sequences = SeqIO.parse(infile, "fasta")
        modified_sequences = []

        for seq_record in sequences:
            seq_record.seq = seq_record.seq.replace("N", "")
            modified_sequences.append(seq_record)
        
        SeqIO.write(modified_sequences, outfile, "fasta")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Remove 'N' characters from sequences in a FASTA file.")
    parser.add_argument("input_file", help="Path to the input FASTA file.")
    parser.add_argument("output_file", help="Path to the output FASTA file.")
    parser.add_argument("fa_dir", help="Path to the output FASTA file.")
    
    args = parser.parse_args()

    remove_N_characters(args.input_file, args.output_file)