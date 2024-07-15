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
doc.margin_x=50
view = GenomeView(chrom, start, end, "+", source)
view.pixel_width = 500
view.pixel_height = 1000

doc.add_view(view)



sample_info=     {
            "Sample": "test_sample",
            "Locus": "DRB1",
            "alleles": ["DRB1*01:01:01:01", "DRB1*01:01:32", "DRB1*01:01:33"],
            "resolution": "4th field"
        }
res_track=ResultTrack(sample_info)
res_track.height = 450
view.margin_x = 50
view.add_track(res_track)
covtrack=BAMCoverageTrack("HLA-C.1.bam", name="bam coverage")
covtrack.margin_y = 50
view.add_track(covtrack)

variant_track = VCFTrack("HG00513.HLA-C.1.phased.vcf.gz", "variants")
print(f"Variant track: {variant_track.var_cnt}")
view.add_track(variant_track)



bam_track_hg002 = SingleEndBAMTrack("HLA-C.1.bam", name="Alignment")
bam_track_hg002.draw_mismatches = True
bam_track_hg002.quick_consensus = True
bam_track_hg002.color_fn = lambda x: "lightgray"
print(bam_track_hg002.mismatch_counts)
view.add_track(bam_track_hg002)

cur_track = BEDTrack('test.bed.gz', name='Gene')
view.add_track(cur_track)

axis_track = Axis()
view.add_track(axis_track)


save(doc, "example.pdf")
