from genomeview import *


def calculate_percent_identity(bam_file):
    bam = pysam.AlignmentFile(bam_file, "rb")
    percent_identity = []
    for read in bam:
        if read.is_unmapped:
            continue
        
        # Calculate the number of identical matches (M) and total alignment length (L)
        num_matches = sum(map(lambda cigar: cigar[1] if cigar[0] == 0 else 0, read.cigartuples))
        total_length = sum(map(lambda cigar: cigar[1] if cigar[0] in (0, 1, 2, 7, 8) else 0, read.cigartuples))
        
        if total_length > 0:
            percent_identity.append((num_matches / total_length) * 100)

            # print(f"{read.query_name}\t{percent_identity:.2f}%")
    
    bam.close()
    return percent_identity[0]

# save(tracks, "example.pdf")
reference_path = "CYP2D6.fasta"
chrom = "CYP2D6"
start = 0
end = 11314
width = 900

source = FastaGenomeSource(reference_path)
doc = Document(width)
doc.margin_x = 50
view = GenomeView(chrom, start, end, "+", source)
view.pixel_width = 500
view.pixel_height = 1000

doc.add_view(view)

sample_info = {
    "Sample": "test_sample",
    "Locus": "DRB1",
    "alleles": ["HLA-DRB3*02:02:01:99000", "HLA-DRB3*02:02:01:99",
                 "HLA-DRB3*02:02:01:09", "HLA-DRB3*02:02:01:09",
                 "HLA-DRB3*02:02:01:09000", "HLA-DRB3*02:02:01:99",
                 "HLA-DRB3*02:02:01:09", "HLA-DRB3*02:02:01:09",
                 "HLA-DRB3*02:02:01:99", "HLA-DRB3*02:02:01:99",
                 "HLA-DRB3*02:02:01:09", "HLA-DRB3*02:02:01:09",
                 "HLA-DRB3*02:02:01:09", "HLA-DRB3*02:02:01:99",
                 "HLA-DRB3*02:02:01:09", "HLA-DRB3*02:02:01:09",
                 ],
    "resolution": "4th field"
}
res_track = ResultTrack(sample_info)

view.add_track(res_track)

covtrack = BAMCoverageTrack("CYP2D6.bam", name="bam coverage")
view.add_track(covtrack)

# variant_track = VCFTrack("HFE.remap2res1.vcf", "variants")
# print(f"Variant track: {variant_track.var_cnt}")
# view.add_track(variant_track)

bam_track_hg002 = SingleEndBAMTrack("CYP2D6.bam", name="Alignment (reference: DRB1*01:01:01:01)", bam_type="normal")
bam_track_hg002.draw_mismatches = True
bam_track_hg002.quick_consensus = True
bam_track_hg002.color_fn = lambda x: "lightgray"
bam_track_hg002.min_indel_size = 20
bam_track_hg002.min_cigar_line_width = 1
print(bam_track_hg002.mismatch_counts)
view.add_track(bam_track_hg002)

# Add another BAM track for a new sample
# bam_identity=calculate_percent_identity("HFE0010101.0.bam")

# show identity in % format

# bam_track_new_sample = SingleEndBAMTrack("HFE0010101.0.bam", name=f"DRB1*01:01:32 (identity: {bam_identity:.2f}%)", bam_type="other")
# bam_track_new_sample.draw_mismatches = True
# bam_track_new_sample.quick_consensus = True
# bam_track_new_sample.color_fn = lambda x: "#588BAF"
# bam_track_new_sample.min_indel_size = 20
# bam_track_new_sample.min_cigar_line_width = 1
# bam_track_new_sample.draw_read_labels = True
# print(bam_track_new_sample.mismatch_counts)
# view.add_track(bam_track_new_sample)

axis_track = Axis()
view.add_track(axis_track)

# cur_track = BEDTrack('test.bed.gz', name='Gene')
# view.add_track(cur_track)



save(doc, "example.pdf")