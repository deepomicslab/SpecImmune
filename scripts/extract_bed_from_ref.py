import os
import subprocess
import shutil
import pysam
import argparse
from Bio import SeqIO

# regions={}
# regions["CYP"]="chr22_KB663609v1_alt chr22_KI270928v1_alt chr22_GL383582v2_alt chr19:40800000-41220000 chr10:93000000-95100000 chr11:14837440-14952231 chr15:51178057-51368601 chr15:74689542-74755536 chr19:15858023-15918077 chr1:46727838-46849413 chr1:47107435-47179735 chr1:59873308-59946773 chr22:42106499-42150865 chr2:38046973-38129902 chr7:963181-1009640 chr7:99600000-99976102 chr7:99235817-99287619 chr20:49483874-49588137"

region={}
region["CYP"]=["chr22_GL383582v2_alt", "chr22_KB663609v1_alt", "chr22_KI270928v1_alt", "chr19:40800000-41220000",
                "chr10:93000000-95100000", "chr11:14837440-14952231", "chr15:51178057-51368601","chr15:74689542-74755536"
               ,"chr19:15858023-15918077","chr1:46727838-46849413","chr1:47107435-47179735","chr1:59873308-59946773",
               "chr22:42106499-42150865","chr2:38046973-38129902","chr7:963181-1009640","chr7:99600000-99976102",
               "chr7:99235817-99287619","chr20:49483874-49588137"]


def align_contigs_to_reference(contigs_file, reference_file, output_bam, threads=1):
    """Align contigs to reference using minimap2 and produce a BAM file."""
    bam_file = output_bam
    # Ensure minimap2 is installed and accessible in PATH
    if not shutil.which("minimap2"):
        raise EnvironmentError("minimap2 is not installed or not in the PATH.")
    
    # Run minimap2 to align the contigs to the reference genome
    cmd = f"minimap2 -t {threads} -a {reference_file} {contigs_file} | samtools view -bS -@ {threads} - | samtools sort -@ {threads} -o {bam_file}"
    subprocess.run(cmd, shell=True, check=True)
    
    # Index the BAM file
    pysam.index(bam_file)
    return bam_file

def extract_reference_regions(bam_file, bed_file, exclude_conventional_chromosomes=False):
    """Extract regions from BAM file based on reference alignment and write to BED file."""
    conventional_chromosomes = {f"chr{i}" for i in range(1, 23)}.union({"chrX", "chrY", "chrM"})
   
    print(region)
    bam = pysam.AlignmentFile(bam_file, "rb")
    aligned_references = set()
    alt_contigs_in_bam=set() 
    with open(bed_file, "w") as bed:
        for read in bam.fetch():
            if not read.is_unmapped:
                chrom = bam.get_reference_name(read.reference_id)
                if exclude_conventional_chromosomes and chrom in conventional_chromosomes:
                    continue
                
                if chrom in aligned_references:
                    continue
                aligned_references.add(chrom)
                alt_contigs_in_bam.add(chrom)
                
                start = 0
                end = bam.get_reference_length(chrom)
                bed.write(f"{chrom}\t{start}\t{end}\n")
                print(f"Writing reference region: {chrom} {start} {end}")

        for default_bed in region["CYP"]:
            bed_items=default_bed.split(":")
            default_contig=bed_items[0]
            print(default_contig)
            if default_contig not in alt_contigs_in_bam:
                if len(bed_items)==1:
                    # write the whole contig to bed
                    bed.write(f"{default_contig}\n")
                    # len_default_contig=bam.get_reference_length(default_contig)
                    print(f"Writing default region: {default_contig}")
                else:
                    start,end=bed_items[1].split("-")
                    bed.write(f"{default_contig}\t{start}\t{end}\n")

                print(f"Writing default region: {default_bed}")

def extract_aligned_regions(bam_file, bed_file, extension_size=100000, exclude_conventional_chromosomes=False):
    """Extract aligned regions from BAM file and write to BED file."""
    conventional_chromosomes = {f"chr{i}" for i in range(1, 23)}.union({"chrX", "chrY", "chrM"})
    
    bam = pysam.AlignmentFile(bam_file, "rb")
    alt_contigs_in_bam=set()
    print(region)
    with open(bed_file, "w") as bed:
        for read in bam.fetch():
            if not read.is_unmapped:
                chrom = bam.get_reference_name(read.reference_id)
                if exclude_conventional_chromosomes and chrom in conventional_chromosomes:
                    continue

                alt_contigs_in_bam.add(chrom)
                start = max(0, read.reference_start - extension_size)
                end = read.reference_end + extension_size
                bed.write(f"{chrom}\t{start}\t{end}\n")
                print(f"Writing aligned region: {chrom} {start} {end}")
        for default_bed in region["CYP"]:
            bed_items=default_bed.split(":")
            default_contig=bed_items[0]
            print(default_contig)
            if default_contig not in alt_contigs_in_bam:
                if len(bed_items)==1:
                    len_default_contig=bam.get_reference_length(default_contig)
                    bed.write(f"{default_contig}\t0\t{len_default_contig}\n")
                    print(f"Writing default region: {default_contig}")
                else:
                    start,end=bed_items[1].split("-")
                    bed.write(f"{default_contig}\t{start}\t{end}\n")

                print(f"Writing default region: {default_bed}")

def merge_bed_regions(input_bed, output_bed):
    """Merge overlapping regions in a BED file using bedtools."""
    # Ensure bedtools is installed and accessible in PATH
    if not shutil.which("bedtools"):
        raise EnvironmentError("bedtools is not installed or not in the PATH.")
    # sort the bed file at first
    cmd = f"bedtools sort -i {input_bed} > {input_bed}.sorted"
    subprocess.run(cmd, shell=True, check=True)
    
    # Run bedtools merge to merge overlapping regions
    cmd = f"bedtools merge -i {input_bed}.sorted > {output_bed}"
    subprocess.run(cmd, shell=True, check=True)

def main(args):
    # Create the working directory if it doesn't exist
    os.makedirs(args.work_dir, exist_ok=True)
    
    intermediate_bed = os.path.join(args.work_dir, "intermediate.bed")
    output_bam = os.path.join(args.work_dir, "aligned.bam")

    # Align contigs to the reference
    bam_file = align_contigs_to_reference(args.contigs_fasta, args.reference_fasta, output_bam, args.threads)
    print(region)
    if args.whole_contigs:
        # Extract whole reference regions to intermediate BED file
        extract_reference_regions(bam_file, intermediate_bed, args.exclude_conventional_chromosomes)
    else:
        # Extract aligned regions to intermediate BED file with extension
        extract_aligned_regions(bam_file, intermediate_bed, args.extension_size, args.exclude_conventional_chromosomes)

    # Merge overlapping regions and write to the final BED file
    if os.path.getsize(intermediate_bed) > 0:
        merge_bed_regions(intermediate_bed, args.work_dir+"/"+args.output_bed)
    else:
        print("No regions to merge; intermediate BED file is empty.")

    # extract bam 

    # Clean up temporary files
    #os.remove(intermediate_bed)
    #os.remove(output_bam)
    #os.remove(output_bam + ".bai")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Align contigs to a reference genome and extract aligned regions.")
    parser.add_argument("contigs_fasta", help="Path to the contigs FASTA file.")
    parser.add_argument("reference_fasta", help="Path to the reference FASTA file.")
    parser.add_argument("output_bed", help="Path to the output BED file.")
    parser.add_argument("-t", "--threads", type=int, default=1, help="Number of threads to use for alignment (default: 1).")
    parser.add_argument("-e", "--exclude_conventional_chromosomes", action="store_true", 
                        help="Exclude conventionally named chromosomes (chr1 to chrM).")
    parser.add_argument("-w", "--whole_contigs", action="store_true", 
                        help="Consider whole reference regions instead of aligned regions.")
    parser.add_argument("-x", "--extension_size", type=int, default=100000, 
                        help="Number of base pairs to extend the aligned regions (default: 100000). Only used if not extracting whole reference regions.")
    parser.add_argument("-d", "--work_dir", default=".", help="Working directory for temporary files (default: current directory).")
    
    args = parser.parse_args()
    main(args)