# from genomeview import *
import sys
sys.path.insert(0, sys.path[0]+'/../../scripts/')
from genomeview import *
import pysam
from Bio import SeqIO


def get_fa_name_and_length(fa_file):
    with open(fa_file, 'r') as f:
        for record in SeqIO.parse(f, "fasta"):
            return record.id, len(record.seq)


ref="HLA-A/HLA-A.fasta"
start = 0
chrom, end = get_fa_name_and_length(ref)
width = 2000

source = FastaGenomeSource(ref)
doc = Document(width)
doc.margin_x = 0
doc.margin_y = 0
view = GenomeView(chrom, start, end, "+", source)
doc.add_view(view)

bam="HLA-A/h0.0.01.bam"
cov_t0=BAMCoverageTrack(bam, name="Bam coverage")
view.add_track(cov_t0)
cov_t0.height = 60

bam_track = SingleEndBAMTrack(bam, name="RNA reads from Haplotype 1", bam_type="normal")
bam_track.draw_mismatches = True
bam_track.quick_consensus = True
bam_track.color_fn = lambda x: "lightgray"
bam_track.min_indel_size = 20
bam_track.min_cigar_line_width = 1
print(bam_track.mismatch_counts)
view.add_track(bam_track)

cur_track = BEDTrack('h0.bed', name="Isoform assembled from haplotype 1 ")
view.add_track(cur_track)


bam="HLA-A/h1.0.01.bam"
cov_t1=BAMCoverageTrack(bam, name="Bam coverage")
cov_t1.height = 60
view.add_track(cov_t1)
bam_track = SingleEndBAMTrack(bam, name="RNA reads from Haplotype 2", bam_type="normal")
bam_track.draw_mismatches = True
bam_track.quick_consensus = True
bam_track.color_fn = lambda x: "lightgray"
bam_track.min_indel_size = 20
bam_track.min_cigar_line_width = 1
print(bam_track.mismatch_counts)
view.add_track(bam_track)

cur_track = BEDTrack('HLA-A/h1.bed', name="Isoform assembled from haplotype 2 ")
view.add_track(cur_track)


axis_track = Axis()
view.add_track(axis_track)


save(doc, 'test.svg')