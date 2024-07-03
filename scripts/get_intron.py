import argparse
from Bio import SeqIO

def parse_gtf(gtf_file):
    intervals = []
    first_isoform_found = False
    transcript_id = None

    with open(gtf_file, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue  # Skip comments
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue  # Skip incomplete lines
            feature_type = parts[2]
            if feature_type != "exon":
                continue  # Skip non-exon features

            # Extract transcript ID from attributes
            attributes = parts[8]
            transcript_id_value = None
            for attribute in attributes.split(';'):
                if 'transcript_id' in attribute:
                    transcript_id_value = attribute.split('"')[1]
                    break

            if not first_isoform_found:
                transcript_id = transcript_id_value
                first_isoform_found = True

            if transcript_id_value != transcript_id:
                continue  # Skip exons not belonging to the first isoform

            seqname = parts[0]
            start = int(parts[3])
            end = int(parts[4])
            intervals.append((start, end))
    return intervals

def merge_intervals(intervals):
    if not intervals:
        return []

    # Sort intervals by starting position
    intervals.sort(key=lambda x: x[0])
    merged = [intervals[0]]
    for current in intervals[1:]:
        last = merged[-1]
        if current[0] <= last[1]:  # Overlapping intervals, merge them
            merged[-1] = (last[0], max(last[1], current[1]))
        else:
            merged.append(current)
    return merged

def find_uncovered_regions(merged_intervals, reference_length):
    uncovered_regions = []
    current_position = 1

    for start, end in merged_intervals:
        if current_position < start:
            uncovered_regions.append((current_position, start - 1))
        current_position = end + 1

    if current_position <= reference_length:
        uncovered_regions.append((current_position, reference_length))

    return uncovered_regions

def write_bed(uncovered_regions, output_file, reference_name):
    with open(output_file, 'w') as file:
        for start, end in uncovered_regions:
            file.write(f"{reference_name}\t{start - 1}\t{end}\n")

def get_reference_length(fasta_file):
    seq_record = SeqIO.read(fasta_file, "fasta")
    return len(seq_record.seq), seq_record.id

def main():
    parser = argparse.ArgumentParser(description="Find uncovered regions in a reference sequence based on GTF file.")
    parser.add_argument("gtf_file", help="Input GTF file")
    parser.add_argument("reference_fasta", help="Input reference FASTA file")
    parser.add_argument("output_file", help="Output BED file for uncovered regions")
    args = parser.parse_args()

    intervals = parse_gtf(args.gtf_file)
    merged_intervals = merge_intervals(intervals)
    reference_length, reference_name = get_reference_length(args.reference_fasta)
    uncovered_regions = find_uncovered_regions(merged_intervals, reference_length)
    write_bed(uncovered_regions, args.output_file, reference_name)

if __name__ == "__main__":
    main()