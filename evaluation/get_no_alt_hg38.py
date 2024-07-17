def is_primary_contig(header):
    """
    Determine if a contig is a primary contig based on its header.
    This function returns True for primary contigs and False for alternative contigs.
    """
    # List of primary contigs in hg38
    primary_contigs = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY', 'chrM']
    # Check if the header matches any of the primary contigs
    contigname = header.split()[0]
    return contigname in primary_contigs

def remove_alternative_contigs(input_fasta, output_fasta):
    """
    Remove alternative contigs from an hg38 FASTA file.
    """
    with open(input_fasta, 'r') as infile, open(output_fasta, 'w') as outfile:
        write_contig = False
        for line in infile:
            if line.startswith('>'):  # FASTA header line
                write_contig = is_primary_contig(line[1:])
            if write_contig:
                outfile.write(line)

# Example usage
input_fasta = '/mnt/d/HLAPro_backup/Nanopore_optimize/data/hg38/GRCh38_full_analysis_set_plus_decoy_hla.fa'
output_fasta = '/mnt/d/HLAPro_backup/Nanopore_optimize/data/hg38/hg38_no_alt.fa'
remove_alternative_contigs(input_fasta, output_fasta)