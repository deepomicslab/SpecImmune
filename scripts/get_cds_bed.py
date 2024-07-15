import pysam

def extract_blocks_to_bed12(sam_file, bed_file):
    samfile = pysam.AlignmentFile(sam_file, "r")
    bed_entries = []
    unique_exons = set()

    for read in samfile.fetch():
        if read.is_unmapped:
            continue
        
        qname = read.query_name  # Query (CDS) sequence name
        rname = samfile.get_reference_name(read.reference_id)  # Reference (gene) sequence name
        blocks = read.get_blocks()  # List of (start, end) tuples for aligned blocks
        strand = '-' if read.is_reverse else '+'
        score = 0  # Placeholder for score, usually 0 or a relevant metric

        # Extract block sizes and start positions
        block_sizes = [end - start for start, end in blocks]
        block_starts = [start - read.reference_start for start, end in blocks]

        # Create a unique identifier for the exon pattern
        exon_pattern = (tuple(block_sizes), tuple(block_starts))

        # Skip if this exon pattern has been seen before for the same reference
        if exon_pattern in unique_exons:
            continue

        # Add the unique exon pattern to the set
        unique_exons.add(exon_pattern)

        # Prepare BED12 format entry
        thick_start = read.reference_start
        thick_end = read.reference_end
        block_count = len(blocks)
        block_sizes_str = ','.join(map(str, block_sizes)) + ','
        block_starts_str = ','.join(map(str, block_starts)) + ','

        bed_entries.append((rname, read.reference_start, read.reference_end, qname, score, strand,
                            thick_start, thick_end, 0, block_count, block_sizes_str, block_starts_str))

    samfile.close()

    with open(bed_file, 'w') as bed:
        for entry in bed_entries:
            bed.write(f"{entry[0]}\t{entry[1]}\t{entry[2]}\t{entry[3]}\t{entry[4]}\t{entry[5]}\t"
                      f"{entry[6]}\t{entry[7]}\t{entry[8]}\t{entry[9]}\t{entry[10]}\t{entry[11]}\n")

def main():
    sam_file = 'cds2gene.2.sam'
    bed_file = 'transcripts.bed'
    extract_blocks_to_bed12(sam_file, bed_file)
    print(f"BED12 file created: {bed_file}")

if __name__ == '__main__':
    main()