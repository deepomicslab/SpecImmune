import subprocess
from pycirclize import Circos
from pycirclize.parser import Gff
import numpy as np

# Load GTF file (parsed using Gff from pycirclize)
gtf_file = "h0.gtf"
gff = Gff(gtf_file)

# Initialize Circos
circos = Circos(sectors={gff.name: gff.range_size})

sector = circos.sectors[0]

# Run samtools depth to get coverage data
bam_file = "h0.bam"  # Replace with your BAM file path
depth_command = ["samtools", "depth", bam_file]

# Run the command and capture the output
result = subprocess.run(depth_command, stdout=subprocess.PIPE, text=True)

# Parse the output of samtools depth
coverage = []
positions = []
for line in result.stdout.strip().split("\n"):
    chrom, pos, cov = line.split("\t")
    positions.append(int(pos))  # The position
    coverage.append(int(cov))   # The coverage at that position

# Create BAM track
vmin = 0
vmax = max(coverage) + 500
bam_track = sector.add_track((85, 150))  # Adjusted range for the BAM track
bam_track.axis(fc="#F3FBF2", ec="none")

# Plot BAM coverage using a filled line plot (similar to IGV)
bam_track.line(positions, coverage, color="blue", lw=0.5)

# Fill the area below the coverage line (to create an IGV-style histogram)
bam_track.fill_between(positions, coverage, np.zeros_like(coverage), color="lightblue")

# Add CDS Track (adjusted to make room for BAM track)
cds_track = sector.add_track((70, 80))  
cds_track.axis(fc="#EEEEEE", ec="none")

# Function to group features by transcript_id manually (since qualifiers are empty)
def group_exons_by_transcript(features):
    transcripts = {}
    for feature in features:
        # Here we can't get transcript_id from qualifiers, so we fallback to a default or manual parsing.
        # If you know the order or context, you might manually assign transcript_id or use `gene_id`.
        transcript_id = "transcript_1"  # Placeholder, adjust based on your knowledge of the data
        if transcript_id not in transcripts:
            transcripts[transcript_id] = []
        transcripts[transcript_id].append(feature)
    return transcripts

# Extract forward strand exons and group them by transcript (for simplicity, assuming one transcript)
f_cds_feats = gff.extract_features("exon", target_strand=1)
f_transcripts = group_exons_by_transcript(f_cds_feats)

# Function to plot exons for each transcript
def plot_transcript_exons(transcripts, track, r_lim, color):
    for transcript_id, exons in transcripts.items():
        # Sort exons by their start position (accessing the start from location)
        exons = sorted(exons, key=lambda x: int(x.location.start))
        exon_starts = [int(exon.location.start) for exon in exons]
        exon_ends = [int(exon.location.end) for exon in exons]
        exon_starts.append(exon_ends[-1])  # Add the last exon end to make the plot continuous
        
        # Draw a line connecting the exons of the same transcript
        # track.line(exon_starts, [r_lim[1]] * len(exon_starts), color=color, lw=1.5)
        
        # Fill between exons to make them look continuous
        track.fill_between(exon_starts, [r_lim[0]] * len(exon_starts), [r_lim[1]/2] * len(exon_starts), color=color, alpha=0.3)
        
        # Add arrows for each exon to show directionality
        for exon in exons:
            track.genomic_features([exon], plotstyle="arrow", r_lim=r_lim, fc=color)

# Plot forward strand exons by transcript
plot_transcript_exons(f_transcripts, cds_track, r_lim=(75, 80), color="salmon")

# Extract reverse strand exons and group them by transcript
r_cds_feats = gff.extract_features("exon", target_strand=-1)
r_transcripts = group_exons_by_transcript(r_cds_feats)
plot_transcript_exons(r_transcripts, cds_track, r_lim=(70, 75), color="skyblue")

# Plot xticks & intervals on inner position
cds_track.xticks_by_interval(
    interval=1000,
    outer=False,
    show_bottom_line=True,
    label_formatter=lambda v: f"{v/1000:.1f} Kb",
    label_orientation="vertical",
    line_kws=dict(ec="grey"),
)

# Plot the figure
fig = circos.plotfig()

# Save the plot as an SVG
fig.savefig("circos_with_bam_linked_exons.svg", format="svg")