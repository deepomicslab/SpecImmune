from genomeview import *

# dataset_paths = ["HLA-C.1.bam"]
# reference = "HLA-C.2.fasta"

# chrom = "HLA-C_ref2"
# start = 0
# end =   3502

# tracks = visualize_data(dataset_paths, chrom, start, end, reference)
# tracks.add_track(BAMCoverageTrack("HLA-C.1.bam", name="pacbio coverage"))



# save(tracks, "example.pdf")
reference_path = "HLA-C.2.fasta"
chrom = "HLA-C_ref2"
start = 0
end = 4304
width = 900

source = FastaGenomeSource(reference_path)
doc = Document(width)
view = GenomeView(chrom, start, end, "+", source)
doc.add_view(view)

view.add_track(BAMCoverageTrack("HLA-C.1.bam", name="bam coverage"))

variant_track = VCFTrack("HG00513.HLA-C.1.phased.vcf.gz", "variants")
print(f"Variant track: {variant_track.var_cnt}")
view.add_track(variant_track)



bam_track_hg002 = SingleEndBAMTrack("HLA-C.1.bam", name="HLA-A hapotype 0")
bam_track_hg002.draw_mismatches = True
bam_track_hg002.quick_consensus = True
bam_track_hg002.color_fn = lambda x: "lightgray"
print(bam_track_hg002.mismatch_counts)
view.add_track(bam_track_hg002)

axis_track = Axis()
view.add_track(axis_track)

save(doc, "example.pdf")
