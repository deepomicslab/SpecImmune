import pysam
import sys

# Define CIGAR operation constants
BAM_CMATCH = 0
BAM_CINS = 1
BAM_CDEL = 2
BAM_CREF_SKIP = 3
BAM_CSOFT_CLIP = 4
BAM_CHARD_CLIP = 5
BAM_CPAD = 6
BAM_CEQUAL = 7
BAM_CDIFF = 8
BAM_CBACK = 9

def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return ''.join(complement[base] for base in reversed(seq))

def extract_clipped_and_inserted_sequences(bam_file, fastq_file, length_threshold, include_reverse_complement):
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        with open(fastq_file, "w") as fastq:
            seen_reads = set()  # To track written read names
            for read in bam:
                if read.query_sequence is None:
                    continue
                
                start_pos = 0  # Initialize start position for each read
                current_position = read.reference_start  # Start at the reference position

                # Check if the read has clipped sequences or insertions
                if read.cigartuples:
                    # Extract clipped sequences
                    if read.cigartuples[0][0] == BAM_CSOFT_CLIP:  # Soft clipping at the start
                        clipped_seq_start = read.query_sequence[:read.cigartuples[0][1]]
                        clipped_qual_start = read.query_qualities[:read.cigartuples[0][1]] if read.query_qualities else [30] * read.cigartuples[0][1]
                        if len(clipped_seq_start) > length_threshold:
                            read_name_start = f"@read_{read.query_name}_start_length_{len(clipped_seq_start)}"
                            if read_name_start not in seen_reads:
                                fastq.write(f"{read_name_start}\n{clipped_seq_start}\n+\n{''.join(chr(q + 33) for q in clipped_qual_start)}\n")
                                seen_reads.add(read_name_start)
                    
                    if read.cigartuples[-1][0] == BAM_CSOFT_CLIP:  # Soft clipping at the end
                        clipped_seq_end = read.query_sequence[-read.cigartuples[-1][1]:]
                        clipped_qual_end = read.query_qualities[-read.cigartuples[-1][1]:] if read.query_qualities else [30] * read.cigartuples[-1][1]
                        if len(clipped_seq_end) > length_threshold:
                            read_name_end = f"@read_{read.query_name}_end_length_{len(clipped_seq_end)}"
                            if read_name_end not in seen_reads:
                                fastq.write(f"{read_name_end}\n{clipped_seq_end}\n+\n{''.join(chr(q + 33) for q in clipped_qual_end)}\n")
                                seen_reads.add(read_name_end)

                    # Extract insertion sequences
                    for cigar_op, length in read.cigartuples:
                        if cigar_op == BAM_CINS:  # Insertion
                            insertion_seq = read.query_sequence[start_pos:start_pos + length]
                            insertion_qual = read.query_qualities[start_pos:start_pos + length] if read.query_qualities else [30] * length
                            if len(insertion_seq) > length_threshold:
                                read_name_insertion = f"@read_{read.query_name}_insertion_{current_position + 1}_length_{len(insertion_seq)}"
                                if read_name_insertion not in seen_reads:
                                    fastq.write(f"{read_name_insertion}\n{insertion_seq}\n+\n{''.join(chr(q + 33) for q in insertion_qual)}\n")
                                    seen_reads.add(read_name_insertion)

                                # Optionally write the reverse complement
                                if include_reverse_complement:
                                    rev_comp_seq = reverse_complement(insertion_seq)
                                    rev_comp_qual = insertion_qual[::-1]  # Reverse the quality scores
                                    read_name_rev_comp = f"@read_{read.query_name}_insertion_{current_position + 1}_rev_comp_length_{len(rev_comp_seq)}"
                                    if read_name_rev_comp not in seen_reads:
                                        fastq.write(f"{read_name_rev_comp}\n{rev_comp_seq}\n+\n{''.join(chr(q + 33) for q in rev_comp_qual)}\n")
                                        seen_reads.add(read_name_rev_comp)

                        # Update positions and start_pos based on the CIGAR operation
                        if cigar_op in (BAM_CMATCH, BAM_CDEL):  # Match or deletion
                            current_position += length
                        if cigar_op in (BAM_CMATCH, BAM_CINS, BAM_CSOFT_CLIP):  # Match, insertion, or soft clip
                            start_pos += length

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python extract_clipped_and_inserted_sequences.py <input_bam> <output_fastq> <length_threshold> <include_reverse_complement>")
        sys.exit(1)

    input_bam = sys.argv[1]  # BAM file path from command line
    output_fastq = sys.argv[2]  # FASTQ file path from command line
    length_threshold = int(sys.argv[3])  # Length threshold from command line
    include_reverse_complement = sys.argv[4].lower() == 'true'  # Boolean flag from command line
    extract_clipped_and_inserted_sequences(input_bam, output_fastq, length_threshold, include_reverse_complement)