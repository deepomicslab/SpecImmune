import pysam

def extract_blocks_to_bed(sam_file, bed_file):
    samfile = pysam.AlignmentFile(sam_file, "r")
    bed_entries = []

    for read in samfile.fetch():
        if read.is_unmapped:
            continue
        
        qname = read.query_name  # Query (CDS) sequence name
        rname = samfile.get_reference_name(read.reference_id)  # Reference (gene) sequence name
        blocks = read.get_blocks()  # List of (start, end) tuples for aligned blocks

        for block in blocks:
            start, end = block
            # BED format is 0-based start, 1-based end
            bed_entries.append((rname, start, end, qname))
    
    samfile.close()

    with open(bed_file, 'w') as bed:
        for entry in bed_entries:
            bed.write(f"{entry[0]}\t{entry[1]}\t{entry[2]}\t{entry[3]}\n")

def main():
    sam_file = 'cds2gene.sam'
    bed_file = 'transcripts.bed'
    extract_blocks_to_bed(sam_file, bed_file)
    print(f"BED file created: {bed_file}")

if __name__ == '__main__':
    main()