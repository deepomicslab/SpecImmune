
import os
import pysam

class LongReadAligner:
    def __init__(self, reference, fastq, alignment_method="minimap2"):
        self.reference = reference
        self.fastq = fastq
        self.alignment_method = alignment_method
        self.indexed = False

    def index_reference(self):
        if self.alignment_method == "minimap2":
            index_file = "reference.mmi"
            os.system(f"minimap2 -d {index_file} {self.reference}")
            self.indexed = True
        elif self.alignment_method == "bwa":
            index_files = [self.reference + ext for ext in [".amb", ".ann", ".bwt", ".pac", ".sa"]]
            if all(os.path.exists(file) for file in index_files):
                self.indexed = True
            else:
                os.system(f"bwa index {self.reference}")
                self.indexed = True
        else:
            raise ValueError("Invalid alignment method specified.")

    def align_reads(self, output_sam):
        if not self.indexed:
            self.index_reference()

        if self.alignment_method == "minimap2":
            index_file = "reference.mmi"
            os.system(f"minimap2 -a {index_file} {self.fastq} -o {output_sam}")
        elif self.alignment_method == "bwa":
            os.system(f"bwa mem -t 8 {self.reference} {self.fastq} > {output_sam}")
        else:
            raise ValueError("Invalid alignment method specified.")

    def sort_bam(self, input_bam, output_bam):
        os.system(f"samtools sort -o {output_bam} {input_bam}")

    def calculate_depth(self, bam_file, output_depth):
        os.system(f"samtools depth -a {bam_file} > {output_depth}")


### a class to accept read type, nanopore, pacbio or pacbio-hifi, add two functions to choose paramter for minimap2 and flyte
class Read_Type:
    def __init__(self, read_type):
        self.read_type = read_type
    
    def get_minimap2_param(self):
        if self.read_type == "nanopore":
            return " -x map-ont "
        elif self.read_type == "pacbio":
            return " -x map-pb "
        elif self.read_type == "pacbio-hifi":
            return " -x map-hifi "
        else:
            raise ValueError("Invalid read type specified.")
    
    def get_flye_param(self):
        if self.read_type == "nanopore":
            return "--nano-raw"
        elif self.read_type == "pacbio":
            return "--pacbio-raw"
        elif self.read_type == "pacbio-hifi":
            return "--pacbio-hifi"
        else:
            raise ValueError("Invalid read type specified.")


def map2db(args, gene, my_db):

    read_type = Read_Type(args["y"])
    minimap_para = read_type.get_minimap2_param()

    outdir = args["o"] + "/" + args["n"]
    sam = outdir + "/" + args["n"] + "." + gene + ".db.sam"
    bam = outdir + "/" + args["n"] + "." + gene + ".db.bam"
    depth_file = outdir + "/" + args["n"] + "." + gene + ".db.depth"
    sort_depth_file = outdir + "/" + args["n"] + "." + gene + ".db.sort.depth.txt"
    # ref={args["f"] }
    # map raw reads to database

    ref = my_db.get_gene_all_alleles(gene)
    # ref="/mnt/d/HLAPro_backup/Nanopore_optimize/SpecHLA/db/HLA/whole/HLA_A.fasta"

    alignDB_order = f"""
    fq={args["r"]}
    
    outdir={args["o"]}/{args["n"]}
    sample={args["n"]}

    fq=$outdir/{gene}.long_read.fq.gz

    seqtk sample $outdir/{gene}.long_read.fq.gz 500 >$outdir/{gene}.long_read.sub.fq
    fq=$outdir/{gene}.long_read.sub.fq

    ref={ref}

    minimap2 -t {args["j"]} {minimap_para} -E 4,2 -p 0.1 -N 100000 -a $ref $fq > {sam}

    # bwa index $ref
    # bwa mem -R '@RG\\tID:foo\\tSM:bar' -a -t {args["j"]} $ref $fq > {sam}

    samtools view -bS -F 0x800  {sam} | samtools sort - >{bam}
    samtools index {bam}
    samtools depth -aa {bam}>{depth_file}
    rm {sam}
    echo alignment done.
    """
    # if the depth_file is not detected 
    # if not os.path.exists(depth_file):
    #     os.system(alignDB_order)
    # else:
    #     print("Depth file is detected.")
    os.system(alignDB_order)

    return bam, depth_file, sort_depth_file

def map2db_blast(args, gene, my_db):

    outdir = args["o"] + "/" + args["n"]
    blast_file = outdir + "/" + args["n"] + "." + gene + ".db.blast.txt"

    alignDB_order = f"""
    outdir={args["o"]}/{args["n"]}
    sample={args["n"]}
  
    fq=$outdir/{gene}.long_read.fq.gz
    seqtk seq -A $fq > $outdir/{gene}.long_read.fasta
    seqtk sample $outdir/{gene}.long_read.fasta 500 >$outdir/{gene}.long_read.sub.fasta
    echo start blastn...
    blastn -query $outdir/{gene}.long_read.sub.fasta -outfmt 7 -out {blast_file} -db {my_db.get_blast_index(gene)}  -num_threads {args["j"]} -max_target_seqs 10000 -max_hsps 1 -task megablast
    echo blast is done.
    """

    os.system(alignDB_order)

    return blast_file


