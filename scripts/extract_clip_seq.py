import pysam
import sys

def extract_clipped_sequences(bam_file, fasta_file):
    # Open the BAM file
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        # Open the FASTA file for writing
        with open(fasta_file, "w") as fasta:
            for read in bam:
                # Check if the read has clipped sequences
                if read.cigartuples:
                    # Extract clipped sequences
                    if read.cigartuples[0][0] == 4:  # Soft clipping at the start
                        clipped_seq = read.query_sequence[:read.cigartuples[0][1]]
                        fasta.write(f">read_{read.query_name}\n{clipped_seq}\n")
                    if read.cigartuples[-1][0] == 4:  # Soft clipping at the end
                        clipped_seq = read.query_sequence[-read.cigartuples[-1][1]:]
                        fasta.write(f">read_{read.query_name}\n{clipped_seq}\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python extract_clipped_sequences.py <input_bam> <output_fasta>")
        sys.exit(1)

    input_bam = sys.argv[1]  # BAM file path from command line
    output_fasta = sys.argv[2]  # FASTA file path from command line
    extract_clipped_sequences(input_bam, output_fasta)