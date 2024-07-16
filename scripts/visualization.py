
import sys
import os
import pysam
import gzip
import argparse
from db_objects import My_db
import csv
import sys
import os
from Bio import SeqIO
from determine_gene import get_focus_gene, get_folder_list
from alignment_modules import Read_Type
from genomeview import *
from PyPDF2 import  PdfMerger
from folder_objects import My_folder


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

def remove_characters(s):
    """
    Remove '*' and ':' characters from the input string.

    Parameters:
    s (str): The input string.

    Returns:
    str: The string with '*' and ':' characters removed.
    """
    s = s.replace('*', '')
    s = s.replace(':', '')
    return s

def get_fa_name_and_length(fa_file):
    with open(fa_file, 'r') as f:
        for record in SeqIO.parse(f, "fasta"):
            return record.id, len(record.seq)

def files_to_pdf(bam, vcf, ref, gene, out_pdf, alleles, alleles_bams=[], tag="hom", allele_idx=0):
    print(f"bam: {bam}, vcf: {vcf}, ref: {ref}, gene: {gene}, out_pdf: {out_pdf}, alleles: {alleles}, alleles_bams: {alleles_bams}")

    sample_info= {
            "Sample": sample,
            "Locus": gene,
            "alleles": alleles,
            "resolution": "4th field"
    }
    # chose top ten alleles
    alleles_bams = alleles_bams[:10]
    alleles = alleles[:10]

    start = 0
    chrom, end = get_fa_name_and_length(ref)
    width = 900

    source = FastaGenomeSource(ref)
    doc = Document(width)
    doc.margin_x = 50
    view = GenomeView(chrom, start, end, "+", source)
    doc.add_view(view)

    res_track=ResultTrack(sample_info)
    res_track.height = 450
    view.add_track(res_track)

    view.add_track(BAMCoverageTrack(bam, name="bam coverage"))

    variant_track = VCFTrack(vcf, "variants")
    print(f"Variant track: {variant_track.var_cnt}")
    view.add_track(variant_track)



    bam_track = SingleEndBAMTrack(bam, name=f"Alignment (reference: constructed {tag} {allele_idx+1})", bam_type="normal")
    bam_track.draw_mismatches = True
    bam_track.quick_consensus = True
    bam_track.color_fn = lambda x: "lightgray"

    bam_track.min_indel_size = 20
    bam_track.min_cigar_line_width = 1
    print(bam_track.mismatch_counts)
    view.add_track(bam_track)

    axis_track = Axis()
    view.add_track(axis_track)

    # remap allele bams
    for allele_idx, allele_bam in enumerate(alleles_bams):
        allele=alleles[allele_idx]
        bam_identity=calculate_percent_identity(allele_bam)
        bam_track = SingleEndBAMTrack(allele_bam, name=f"{allele} (identity: {bam_identity:.2f}%)", bam_type="parwise")
        bam_track.draw_mismatches = True
        bam_track.quick_consensus = True
        bam_track.color_fn = lambda x: "#588BAF"

        bam_track.min_indel_size = 20
        bam_track.min_cigar_line_width = 1
        print(bam_track.mismatch_counts)
        view.add_track(bam_track)

    save(doc, out_pdf)


def call_variants(bam, ref, ovcf):
            # longshot -F -c 2 -C 100000 -P {args["strand_bias_pvalue_cutoff"]} -r {interval_dict[gene]} --bam {bam} --ref {hla_ref} --out {parameter.outdir}/{parameter.sample}.{gene}.longshot.vcf 
    cmd=f"""
    longshot -F -c 2 -C 100000 -P 0.01 --bam {bam} --ref {ref} --out {ovcf} 
    """
    os.system(cmd)

def read_hla_file(filename, some_dict):
    with open(filename, 'r') as file:
        next(file)
        next(file)
        for idx, row in enumerate(file):
            # print(row.strip())
            items=row.strip().split("\t")
            if len(items)<3:
                # remove gene_tag from some_dict
                gene_tag = items[0]
                some_dict[gene_tag].append("-")
                continue
            elif items[2] == "":
                continue
            alleles_str=items[2].split(";")
            print("alleles_str", alleles_str)
            for iti, it in enumerate(alleles_str):
                if iti>0:
                    continue
                allele=it.split(",")[0]
                # gene_tag=allele.split("*")[0]
                gene_tag = items[0]
                print("result :", gene_tag, allele) 
                some_dict[gene_tag].append(allele)

def parse_full_allele(filename, full_allele_dict):
    with open(filename, 'r') as file:
        next(file)
        next(file)
        for idx, row in enumerate(file):
            print(row.strip())
            items=row.strip().split("\t")
            if len(items)<3:
                # remove gene_tag from some_dict
                gene_tag = items[0]
                continue
            elif items[2] == "":
                continue
        
            alleles_str=items[2].split(";")
            print("alleles_str", alleles_str)
            for iti, it in enumerate(alleles_str):
                allele=it.split(",")[0]
                # gene_tag=allele.split("*")[0]
                gene_tag = items[0]
                gene_idx=int(items[1])-1
                print("result :", gene_tag, allele) 
                full_allele_dict[gene_tag][gene_idx].append(allele)

def generate_pdf():
    for gene, alleles in step1_res_dict.items():
        gene_dir=f"{remap_allele_dir}/{gene}"
        if len(alleles)==0:
            continue
        else:
            if '-' == alleles[0]:
                continue
    # for step1 hom
        if alleles[0] == alleles[1]:
            # for step2 hom
            if len(step2_res_dict[gene]) == 0:
                continue
            else:
                if "-" == step2_res_dict[gene][0]:
                    continue
            if step2_res_dict[gene][0] == step2_res_dict[gene][1]:
                ref=f"{my_folder.sequence_dir}/{gene_class}.allele.1.{gene}.fasta"
                bam=f"{my_folder.for_viz_dir}/{gene}.remap2res.bam"
                vcf=f"{my_folder.for_viz_dir}/{gene}.remap2res.vcf"
                out_pdf=f"{my_folder.visualization_dir}/{gene}.remap2res.pdf"
                allele_remap_bams=[f"{gene_dir}/{remove_characters(allele)}.{allele_idx}.bam" for allele in full_allele_dict[gene][0]]
                if os.path.exists(ref):
                    files_to_pdf(bam, vcf, ref, gene, out_pdf, full_allele_dict[gene][0], allele_remap_bams, "hom", 0)
                
            else:
                # for step2 het
                ref1=f"{my_folder.sequence_dir}/{gene_class}.allele.1.{gene}.fasta"
                ref2=f"{my_folder.sequence_dir}/{gene_class}.allele.2.{gene}.fasta"
                bam1=f"{my_folder.for_viz_dir}/{gene}.remap2res1.bam"
                bam2=f"{my_folder.for_viz_dir}/{gene}.remap2res2.bam"
                vcf1=f"{my_folder.for_viz_dir}/{gene}.remap2res1.vcf"
                vcf2=f"{my_folder.for_viz_dir}/{gene}.remap2res2.vcf"
                out_pdf1=f"{my_folder.visualization_dir}/{gene}.remap2res1.pdf"
                out_pdf2=f"{my_folder.visualization_dir}/{gene}.remap2res2.pdf"
                # ref1=f"{outdir}/{sample}/{gene_class}.allele.1.{gene}.fasta"
                # ref2=f"{outdir}/{sample}/{gene_class}.allele.2.{gene}.fasta"
                # bam1=f"{outdir}/{sample}/{gene}.remap2res1.bam"
                # bam2=f"{outdir}/{sample}/{gene}.remap2res2.bam"
                # vcf1=f"{outdir}/{sample}/{gene}.remap2res1.vcf"
                # vcf2=f"{outdir}/{sample}/{gene}.remap2res2.vcf"
                # out_pdf1=f"{outdir}/{sample}/{gene}.remap2res1.pdf"
                # out_pdf2=f"{outdir}/{sample}/{gene}.remap2res2.pdf"
                if os.path.exists(ref1):
                    allele_remap_bams1=[f"{gene_dir}/{remove_characters(allele)}.{0}.bam" for allele in full_allele_dict[gene][0]]
                    files_to_pdf(bam1, vcf1, ref1, gene, out_pdf1, full_allele_dict[gene][0], allele_remap_bams1, "het", 0)
                if os.path.exists(ref2):
                    allele_remap_bams2=[f"{gene_dir}/{remove_characters(allele)}.{1}.bam" for allele in full_allele_dict[gene][1]]
                    files_to_pdf(bam2, vcf2, ref2, gene, out_pdf2, full_allele_dict[gene][1], allele_remap_bams2, "het", 1)

        else:
            # for step1 het
            for allele_idx, allele in enumerate(alleles):
                # for step2 hom
                if len(step2_res_dict[gene]) == 0:
                    continue
                else:
                    if "-" == step2_res_dict[gene][0]:
                        continue
                if step2_res_dict[gene][0] == step2_res_dict[gene][1]:
                    ref=f"{my_folder.sequence_dir}/{gene_class}.allele.1.{gene}.fasta"
                    bam=f"{my_folder.for_viz_dir}/{gene}.remap2res{allele_idx+1}.bam"
                    vcf=f"{my_folder.for_viz_dir}/{gene}.remap2res{allele_idx+1}.vcf"
                    out_pdf=f"{my_folder.visualization_dir}/{gene}.remap2res{allele_idx+1}.pdf"

                    if os.path.exists(ref):
                        allele_remap_bams=[f"{gene_dir}/{remove_characters(allele)}.{allele_idx}.bam" for allele in full_allele_dict[gene][allele_idx]]
                        files_to_pdf(bam, vcf, ref, gene, out_pdf, full_allele_dict[gene][allele_idx], allele_remap_bams, "hom", allele_idx)
                else:
                    # for step2 het
                    ref1=f"{my_folder.sequence_dir}/{gene_class}.allele.1.{gene}.fasta"
                    ref2=f"{my_folder.sequence_dir}/{gene_class}.allele.2.{gene}.fasta"
                    bam1=f"{my_folder.for_viz_dir}/{gene}.remap2res1.bam"
                    bam2=f"{my_folder.for_viz_dir}/{gene}.remap2res2.bam"
                    vcf1=f"{my_folder.for_viz_dir}/{gene}.remap2res1.vcf"
                    vcf2=f"{my_folder.for_viz_dir}/{gene}.remap2res2.vcf"
                    out_pdf1=f"{my_folder.visualization_dir}/{gene}.remap2res1.pdf"
                    out_pdf2=f"{my_folder.visualization_dir}/{gene}.remap2res2.pdf"
                
                    if os.path.exists(ref1):
                        allele_remap_bams1=[f"{gene_dir}/{remove_characters(allele)}.{0}.bam" for allele in full_allele_dict[gene][0]]
                        print("allele_remap_bams1", allele_remap_bams1)
                        files_to_pdf(bam1, vcf1, ref1, gene, out_pdf1, full_allele_dict[gene][0], allele_remap_bams1, "het", 0)
                    if os.path.exists(ref2):
                        allele_remap_bams2=[f"{gene_dir}/{remove_characters(allele)}.{1}.bam" for allele in full_allele_dict[gene][1]]
                        print("allele_remap_bams2", allele_remap_bams2)
                        files_to_pdf(bam2, vcf2, ref2, gene, out_pdf2, full_allele_dict[gene][1], allele_remap_bams2, "het", 1)
                    break

def merge_pdfs(pdf_list, output_path):
    pdf_merger = PdfMerger()

    for pdf in pdf_list:
        try:
            if not os.path.isfile(pdf):
                print(f"File not found: {pdf}")
                continue
            pdf_merger.append(pdf)
        except PyPDF2.errors.PdfReadError as e:
            print(f"Error reading {pdf}: {e}")
        except Exception as e:
            print(f"An unexpected error occurred with {pdf}: {e}")

    with open(output_path, 'wb') as output_file:
        pdf_merger.write(output_file)

def main():
    # parse step one res
    read_hla_file(step1_result, step1_res_dict)
    read_hla_file(step2_result, step2_res_dict)
    parse_full_allele(step2_result, full_allele_dict)
    generate_pdf()
    # concatenate all the pdfs to one pdf
    # find *pdf files in the outdir
    pdf_files = [f"{my_folder.visualization_dir}/{f}" for f in os.listdir(f"{my_folder.visualization_dir}") if f.endswith(".pdf")]

    merge_pdfs(pdf_files, f"{outdir}/{sample}/{sample}.pdf")


if __name__ == "__main__":   
    if len(sys.argv) <2:
        print("Usage: python script.py <filename>")
        sys.exit(1)
    
    sample = sys.argv[1]
    gene_class = sys.argv[2]
    outdir = sys.argv[3]
    data_type = sys.argv[4]
    seq_tech = sys.argv[5]
    RNA_type = sys.argv[6]
    threads = sys.argv[7]
    db_ref = sys.argv[8]

    read_type = Read_Type(seq_tech, data_type, RNA_type)
    my_folder = My_folder({"o": outdir, "n":sample})

    step1_result = f"{my_folder.sample_prefix}.{gene_class}.type.result.txt"
    step2_result = f"{my_folder.sample_prefix}.{gene_class}.final.type.result.txt"   #f"{my_folder.outdir}/hlala.like.results.txt"

    remap_allele_dir=my_folder.for_viz_dir + "/remap_allele"

    minimap_para = read_type.get_minimap2_param()
    step1_res_dict = {}
    step2_res_dict = {}
    full_allele_dict = {}
    db_folder=os.path.dirname(db_ref)
    gene_list = get_folder_list(db_folder)
    for gene in gene_list:
        step1_res_dict[gene] = []
        step2_res_dict[gene] = []
        full_allele_dict[gene] = [[],[]]
    main()
        



