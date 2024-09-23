import os
import subprocess
from Bio import SeqIO  # Import SeqIO from Biopython to handle FASTA files
from pycirclize import Circos
from pycirclize.parser import Gff
import numpy as np

# List of genes and their respective directories
genes = ["HLA-A", "HLA-B"]
base_dir = "./"  # Adjust this path

# Dictionary to store the sectors (gene names with their range sizes)
sectors = {}

# Iterate over each gene to gather contig length from its respective FASTA file
for gene in genes:
    gene_dir = os.path.join(base_dir, gene)  # Path to the gene's directory
    fasta_file = os.path.join(gene_dir, f"{gene}.fasta")  # Path to the gene's FASTA file

    # Check if the FASTA file exists
    if not os.path.isfile(fasta_file):
        raise FileNotFoundError(f"FASTA file for {gene} not found in {gene_dir}")

    # Parse the FASTA file and extract the contig length (assuming one contig per gene)
    with open(fasta_file, "r") as fasta_handle:
        for record in SeqIO.parse(fasta_handle, "fasta"):
            contig_length = len(record.seq)
            sectors[gene] = contig_length  # Store the contig length for this gene
            break  # Assuming only one contig per gene FASTA file

# Now, initialize Circos with the sectors (genes with their contig lengths)
circos = Circos(sectors=sectors, space=10)

# Function to group exons by transcript (used later)
def group_exons_by_transcript(features):
    transcripts = {}
    for feature in features:
        transcript_id = "transcript_1"  # Placeholder; replace based on your knowledge of the data
        if transcript_id not in transcripts:
            transcripts[transcript_id] = []
        transcripts[transcript_id].append(feature)
    return transcripts

# Function to plot exons for each transcript (used later)
def plot_transcript_exons(transcripts, track, r_lim, color, is_forward=True):
    for transcript_id, exons in transcripts.items():
        # Sort exons by start position
        exons = sorted(exons, key=lambda x: int(x.location.start))
        exon_starts = [int(exon.location.start) for exon in exons]
        exon_ends = [int(exon.location.end) for exon in exons]
        exon_starts.append(exon_ends[-1])  # Add the last exon end to make the plot continuous

        # Fill between the exons with different heights depending on strand direction
        if is_forward:
            print(r_lim,[r_lim[0]] * len(exon_starts), [r_lim[1] / 2] * len(exon_starts) )
            track.fill_between(
                exon_starts, 
                [r_lim[0]] * len(exon_starts), 
                [r_lim[1] / 2] * len(exon_starts),  
                color=color, alpha=0.3
            )
        else:
            print(r_lim,[r_lim[0]] * len(exon_starts), [r_lim[1]] * len(exon_starts) )
            track.fill_between(
                x=exon_starts, 
                y1=[60] * len(exon_starts), 
                y2=[75] * len(exon_starts),  
                color=color, alpha=0.3
            )

        # Plot exons as arrows
        for exon in exons:
            track.genomic_features([exon], plotstyle="arrow", r_lim=r_lim, fc=color)

# Iterate over each gene to plot the data
for idx, gene in enumerate(genes):  # Use enumerate to get both index and gene name
    gene_dir = os.path.join(base_dir, gene)  # Path to the gene's directory
    gtf_file = os.path.join(gene_dir, "h0.gff")  # Path to the GFF file
    bam_file = os.path.join(gene_dir, "h0.bam")  # Path to the BAM file

    # Load GFF file
    gff = Gff(gtf_file)

    # Get the sector by index (use idx to access the correct sector)
    sector = circos.sectors[idx]  # Access sectors by index

    # Add a track for the sector label (gene name) using text
    label_track = sector.add_track((95, 100))  # Add the label at an outer position (adjust range as needed)
    label_track.text(gene, size=12, color="black", orientation="vertical", va="center")  # Adjust size, color, and orientation

    # Run samtools depth to get coverage data from the BAM file
    depth_command = ["samtools", "depth", bam_file]
    result = subprocess.run(depth_command, stdout=subprocess.PIPE, text=True)

    # Parse samtools depth output
    coverage = []
    positions = []
    for line in result.stdout.strip().split("\n"):
        chrom, pos, cov = line.split("\t")
        positions.append(int(pos))
        coverage.append(int(cov))

    # Plot BAM coverage
    vmin = 0
    vmax = max(coverage) + 500
    bam_track = sector.add_track((85, 100))
    bam_track.axis(fc="#F3FBF2", ec="lightblue", lw=0.5)  # Add border with black color and width of 1.5
    bam_track.line(positions, coverage, color="blue", lw=0.5)
    bam_track.fill_between(positions, coverage, np.zeros_like(coverage), color="lightblue")

    # Plot CDS (exon) track
    cds_track = sector.add_track((70, 80))
    cds_track.axis(fc="#EEEEEE", ec="lightblue", lw=1)  # Add border with black color and width of 1.5

    # Extract forward strand exons and group them by transcript
    f_cds_feats = gff.extract_features("exon", target_strand=1)
    f_transcripts = group_exons_by_transcript(f_cds_feats)
    plot_transcript_exons(f_transcripts, cds_track, r_lim=(75, 80), color="salmon", is_forward=True)

    # Extract reverse strand exons and group them by transcript
    r_cds_feats = gff.extract_features("exon", target_strand=-1)
    r_transcripts = group_exons_by_transcript(r_cds_feats)
    print(r_transcripts)
    plot_transcript_exons(r_transcripts, cds_track, r_lim=(70, 75), color="skyblue", is_forward=False)

    # Plot xticks & intervals on inner position for this sector
    cds_track.xticks_by_interval(
        interval=1000,
        outer=False,
        show_bottom_line=True,
        label_formatter=lambda v: f"{v / 1000:.1f} Kb",
        label_orientation="vertical",
        line_kws=dict(ec="grey"),
    )

# Finally, plot the Circos figure
fig = circos.plotfig()

# Save the plot as an SVG
fig.savefig("circos_with_multiple_genes.png", format="png")