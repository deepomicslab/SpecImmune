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



<<<<<<< HEAD
fa_allele_dir="c_HG00420/alleles/fa_alleles/"
ma_allele_dir="c_HG00420/alleles/ma_alleles/"

bam_dir="c_HG00420/remap/"
vcf_dir="c_HG00420/remap/"
outdir="c_HG00420/viz/"
=======
fa_allele_dir="c_NA19828/alleles/fa_alleles/"
ma_allele_dir="c_NA19828/alleles/ma_alleles/"

bam_dir="c_NA19828/remap/"
vcf_dir="c_NA19828/remap/"
outdir="c_NA19828/viz/"
>>>>>>> 332fc61c7a58d91ab791aa4db2429f96e21f3c2c
if not os.path.exists(outdir):
    os.makedirs(outdir)

gene_list = [ 'HLA-A', 'HLA-B', 'HLA-C', 'HLA-DMA', 'HLA-DMB', 'HLA-DOA', 'HLA-DOB', 'HLA-DPA1', 'HLA-DPB1', 'HLA-DPB2', 'HLA-DQA1', 'HLA-DQB1','HLA-DQB2', 'HLA-DRA', 'HLA-DRB1', 'HLA-DRB3', 'HLA-DRB4', 'HLA-DRB5', 'HLA-E', 'HLA-F', 'HLA-G', 'HLA-H', 'HLA-J', 'HLA-K', 'HLA-L', 'HLA-P', 'HLA-V', 'HLA-Y', 'HLA-DQA2', 'HLA-DPA2', 'HLA-N', 'HLA-S', 'HLA-T', 'HLA-U', 'HLA-W', 'MICA', 'MICB', 'TAP1', 'TAP2', 'HFE' ]


start = 0
width = 1000
for gene in gene_list:
    fa_allele1=f"{fa_allele_dir}/HLA.allele.1.{gene}.fasta"
    c_map_fa_allele1_vcf=f"{bam_dir}/c_{gene}_fa1.vcf"
    c_map_fa_allele1_bam=f"{vcf_dir}/c_{gene}_fa1.bam"
    c_reads_map_fa_allele1_bam=f"{vcf_dir}/reads_{gene}_fa1.bam"

    # fa_allele1 or c_map_fa_allele1_vcf or c_map_fa_allele1_bam or c_reads_map_fa_allele1_bam not found, skip
    if not os.path.exists(fa_allele1) or not os.path.exists(c_map_fa_allele1_vcf) or not os.path.exists(c_map_fa_allele1_bam) or not os.path.exists(c_reads_map_fa_allele1_bam):
        continue

    print(f"c_map_fa_allele1_vcf: {c_map_fa_allele1_vcf}")
    print(f"c_map_fa_allele1_bam: {c_map_fa_allele1_bam}")
    chrom, end = get_fa_name_and_length(fa_allele1)
    source = FastaGenomeSource(fa_allele1)
    doc = Document(width)
    doc.margin_x = 0
    doc.margin_y = 0
    view = GenomeView(chrom, start, end, "+", source)
    doc.add_view(view)


    variant_track = VCFTrack(c_map_fa_allele1_vcf, "variants")
    print(f"Variant track: {variant_track.var_cnt}")
    view.add_track(variant_track)


    cov=BAMCoverageTrack(f"{c_reads_map_fa_allele1_bam}", name="Bam coverage")
    view.add_track(cov)
    # cov.height = 70

    bam_track = SingleEndBAMTrack(f"{c_reads_map_fa_allele1_bam}", name="Child reads align to father allele1", bam_type="normal")
    bam_track.draw_mismatches = True
    bam_track.quick_consensus = True
    bam_track.color_fn = lambda x: "lightgray"
    bam_track.min_indel_size = 20
    bam_track.min_cigar_line_width = 1
    # bam_track.height = 50
    view.add_track(bam_track)


    cov=BAMCoverageTrack(f"{c_map_fa_allele1_bam}", name="Bam coverage")
    view.add_track(cov)
    # cov.height = 50

    bam_track = SingleEndBAMTrack(f"{c_map_fa_allele1_bam}", name="Child haplotype align to father allele1", bam_type="parwise")
    bam_track.draw_mismatches = True
    bam_track.quick_consensus = True
    bam_track.color_fn = lambda x: "lightgray"
    bam_track.min_indel_size = 20
    bam_track.min_cigar_line_width = 1
<<<<<<< HEAD
    # bam_track.height = 70
=======
    bam_track.height = 70
>>>>>>> 332fc61c7a58d91ab791aa4db2429f96e21f3c2c
    view.add_track(bam_track)

    fa_allele2=f"{fa_allele_dir}/HLA.allele.2.{gene}.fasta"
    c_map_fa_allele2_vcf=f"{bam_dir}/c_{gene}_fa2.vcf"
    c_map_fa_allele2_bam=f"{vcf_dir}/c_{gene}_fa2.bam"
    c_reads_map_fa_allele2_bam=f"{vcf_dir}/reads_{gene}_fa2.bam"
    if not os.path.exists(fa_allele2) or not os.path.exists(c_map_fa_allele2_vcf) or not os.path.exists(c_map_fa_allele2_bam) or not os.path.exists(c_reads_map_fa_allele2_bam):
        continue
    print(f"c_map_fa_allele2_vcf: {c_map_fa_allele2_vcf}")
    print(f"c_map_fa_allele2_bam: {c_map_fa_allele2_bam}")
    chrom, end = get_fa_name_and_length(fa_allele2)
    source = FastaGenomeSource(fa_allele2)
    view = GenomeView(chrom, start, end, "+", source)
    doc.add_view(view)


    variant_track = VCFTrack(c_map_fa_allele2_vcf, "variants")
    print(f"Variant track: {variant_track.var_cnt}")
    view.add_track(variant_track)


    cov=BAMCoverageTrack(f"{c_reads_map_fa_allele2_bam}", name="Bam coverage")
    view.add_track(cov)
<<<<<<< HEAD
    # cov.height = 70
=======
    cov.height = 70
>>>>>>> 332fc61c7a58d91ab791aa4db2429f96e21f3c2c

    bam_track = SingleEndBAMTrack(f"{c_reads_map_fa_allele2_bam}", name="Child reads align to father allele2", bam_type="normal")
    bam_track.draw_mismatches = True
    bam_track.quick_consensus = True
    bam_track.color_fn = lambda x: "lightgray"
    bam_track.min_indel_size = 20
    bam_track.min_cigar_line_width = 1
    # bam_track.height = 50
    view.add_track(bam_track)


    cov=BAMCoverageTrack(f"{c_map_fa_allele2_bam}", name="Bam coverage")
    view.add_track(cov)
    # cov.height = 50

    bam_track = SingleEndBAMTrack(f"{c_map_fa_allele2_bam}", name="Child haplotype align to father allele2", bam_type="parwise")
    bam_track.draw_mismatches = True
    bam_track.quick_consensus = True
    bam_track.color_fn = lambda x: "lightgray"
    bam_track.min_indel_size = 20
    bam_track.min_cigar_line_width = 1
<<<<<<< HEAD
    # bam_track.height = 70
=======
    bam_track.height = 70
>>>>>>> 332fc61c7a58d91ab791aa4db2429f96e21f3c2c
    view.add_track(bam_track)



    ma_allele1=f"{ma_allele_dir}/HLA.allele.1.{gene}.fasta"
    c_map_ma_allele1_vcf=f"{bam_dir}/c_{gene}_ma1.vcf"
    c_map_ma_allele1_bam=f"{vcf_dir}/c_{gene}_ma1.bam"
    c_reads_map_ma_allele1_bam=f"{vcf_dir}/reads_{gene}_ma1.bam"
    if not os.path.exists(ma_allele1) or not os.path.exists(c_map_ma_allele1_vcf) or not os.path.exists(c_map_ma_allele1_bam) or not os.path.exists(c_reads_map_ma_allele1_bam):
        continue
    print(f"c_map_ma_allele1_vcf: {c_map_ma_allele1_vcf}")
    print(f"c_map_ma_allele1_bam: {c_map_ma_allele1_bam}")
    chrom, end = get_fa_name_and_length(ma_allele1)
    source = FastaGenomeSource(ma_allele1)
    view = GenomeView(chrom, start, end, "+", source)
    doc.add_view(view)


    variant_track = VCFTrack(c_map_ma_allele1_vcf, "variants")
    print(f"Variant track: {variant_track.var_cnt}")
    view.add_track(variant_track)


    cov=BAMCoverageTrack(f"{c_reads_map_ma_allele1_bam}", name="Bam coverage")
    view.add_track(cov)
    # cov.height = 70

    bam_track = SingleEndBAMTrack(f"{c_reads_map_ma_allele1_bam}", name="Child reads align to mother allele1", bam_type="normal")
    bam_track.draw_mismatches = True
    bam_track.quick_consensus = True
    bam_track.color_fn = lambda x: "lightgray"
    bam_track.min_indel_size = 20
    bam_track.min_cigar_line_width = 1
    # bam_track.height = 50
    view.add_track(bam_track)


    cov=BAMCoverageTrack(f"{c_map_ma_allele1_bam}", name="Bam coverage")
    view.add_track(cov)
    # cov.height = 50

    bam_track = SingleEndBAMTrack(f"{c_map_ma_allele1_bam}", name="Child haplotype align to mother allele1", bam_type="parwise")
    bam_track.draw_mismatches = True
    bam_track.quick_consensus = True
    bam_track.color_fn = lambda x: "lightgray"
    bam_track.min_indel_size = 20
    bam_track.min_cigar_line_width = 1
<<<<<<< HEAD
    # bam_track.height = 70
=======
    bam_track.height = 70
>>>>>>> 332fc61c7a58d91ab791aa4db2429f96e21f3c2c
    view.add_track(bam_track)

    ma_allele2=f"{ma_allele_dir}/HLA.allele.2.{gene}.fasta"
    c_map_ma_allele2_vcf=f"{bam_dir}/c_{gene}_ma2.vcf"
    c_map_ma_allele2_bam=f"{vcf_dir}/c_{gene}_ma2.bam"
    c_reads_map_ma_allele2_bam=f"{vcf_dir}/reads_{gene}_ma2.bam"
    if not os.path.exists(ma_allele2) or not os.path.exists(c_map_ma_allele2_vcf) or not os.path.exists(c_map_ma_allele2_bam) or not os.path.exists(c_reads_map_ma_allele2_bam):
        continue
    print(f"c_map_ma_allele2_vcf: {c_map_ma_allele2_vcf}")
    print(f"c_map_ma_allele2_bam: {c_map_ma_allele2_bam}")
    chrom, end = get_fa_name_and_length(ma_allele2)
    source = FastaGenomeSource(ma_allele2)
    view = GenomeView(chrom, start, end, "+", source)
    doc.add_view(view)


    variant_track = VCFTrack(c_map_ma_allele2_vcf, "variants")
    print(f"Variant track: {variant_track.var_cnt}")
    view.add_track(variant_track)

    
    cov=BAMCoverageTrack(f"{c_reads_map_ma_allele2_bam}", name="Bam coverage")
    view.add_track(cov)
<<<<<<< HEAD
    # cov.height = 70
=======
    cov.height = 70
>>>>>>> 332fc61c7a58d91ab791aa4db2429f96e21f3c2c

    bam_track = SingleEndBAMTrack(f"{c_reads_map_ma_allele2_bam}", name="Child reads align to mother allele2", bam_type="normal")
    bam_track.draw_mismatches = True
    bam_track.quick_consensus = True
    bam_track.color_fn = lambda x: "lightgray"
    bam_track.min_indel_size = 20
    bam_track.min_cigar_line_width = 1
    # bam_track.height = 50
    view.add_track(bam_track)


    cov=BAMCoverageTrack(f"{c_map_ma_allele2_bam}", name="Bam coverage")
    view.add_track(cov)
    # cov.height = 50

    bam_track = SingleEndBAMTrack(f"{c_map_ma_allele2_bam}", name="Child haplotype align to mother allele2", bam_type="parwise")
    bam_track.draw_mismatches = True
    bam_track.quick_consensus = True
    bam_track.color_fn = lambda x: "lightgray"
    bam_track.min_indel_size = 20
    bam_track.min_cigar_line_width = 1
<<<<<<< HEAD
    # bam_track.height = 70
=======
    bam_track.height = 70
>>>>>>> 332fc61c7a58d91ab791aa4db2429f96e21f3c2c
    view.add_track(bam_track)












    save(doc, f'{outdir}/{gene}.pdf')





