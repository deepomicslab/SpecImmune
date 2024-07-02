
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
    def __init__(self, seq_tech, read_type, RNA_type):
        self.seq_tech=seq_tech
        self.read_type = read_type
        self.RNA_type = RNA_type
    
    def get_minimap2_param(self):
        if self.seq_tech in ["wgs", "amplicon"]:
            if self.read_type == "nanopore":
                return " -x map-ont "
            elif self.read_type == "pacbio":
                return " -x map-pb "
            elif self.read_type == "pacbio-hifi":
                return " -x map-hifi "
            else:
                raise ValueError("Invalid read type specified.")
        # for splice alignment
        elif self.seq_tech in ["rna"]:
            if self.RNA_type == "traditional":
                return " -x splice:hq -uf "
            elif self.RNA_type == "2D":
                return " -x splice "
            elif self.RNA_type == "Direct":
                return " -x splice -uf -k14 "
            elif self.RNA_type == "SIRV":
                return " -x splice --splice-flank=no "
            else:
                raise ValueError("Invalid RNA type specified.")
            
    # for cds alignment
    def get_bwa_param(self):
        if self.seq_tech in ["rna"]:
            if self.RNA_type == "traditional":
                return " -x pacbio "
            elif self.RNA_type == "2D":
                return " -x ont2d "
            elif self.RNA_type == "Direct":
                return " -x pacbio "
            elif self.RNA_type == "SIRV":
                return " -x pacbio "
            else:
                raise ValueError("Invalid RNA type specified.")

    
    def get_flye_param(self):
        if self.read_type == "nanopore":
            return "--nano-raw"
        elif self.read_type == "pacbio":
            return "--pacbio-raw"
        elif self.read_type == "pacbio-hifi":
            return "--pacbio-hifi"
        else:
            raise ValueError("Invalid read type specified.")

def cout_read_num(fastq): 
    ## the input fastq is gziped, count the reads num in it, and return it to the variable
    f = open(fastq, "rb")
    count = 0
    for line in f:
        count += 1
    f.close()
    return count/4

def subsample_fastq(gene, args, read_num):
    outdir = args["o"] + "/" + args["n"]
    fastq = outdir + "/" + gene + ".long_read.fq.gz"
    sub_fastq = outdir + "/" + gene + ".long_read.sub.fq"
    read_count = cout_read_num(fastq)
    if read_count <= read_num:
        os.system(f"zcat {fastq} > {sub_fastq}")
    else:
        
        
        os.system(f"seqtk sample {fastq} {read_num} > {sub_fastq}")
        ## gzip the fq
        # os.system(f"gzip -f {sub_fastq}")
    return sub_fastq





def map2db(args, gene, my_db, read_num=500):
    # map binned reads to all alleles of each locus
    read_type = Read_Type(args["seq_tech"], args["y"], args["RNA_type"])
    if args["seq_tech"] == "rna":
        bwa_para = read_type.get_bwa_param()
    else:
        minimap_para = read_type.get_minimap2_param()
    # minimap_para = read_type.get_minimap2_param()

    outdir = args["o"] + "/" + args["n"]
    sam = outdir + "/" + args["n"] + "." + gene + ".db.sam"
    bam = outdir + "/" + args["n"] + "." + gene + ".db.bam"
    depth_file = outdir + "/" + args["n"] + "." + gene + ".db.depth"
    sort_depth_file = outdir + "/" + args["n"] + "." + gene + ".db.sort.depth.txt"
    # ref={args["f"] }
    # map raw reads to database

    ref = my_db.get_gene_all_alleles(gene)
    sub_fastq = subsample_fastq(gene, args, read_num)

    if args["seq_tech"] == "rna":
        alignDB_order = f"""
        outdir={args["o"]}/{args["n"]}
        sample={args["n"]}
        ref={ref}
        bwa mem {bwa_para} -t {args["j"]} $ref {sub_fastq} | samtools view -bS -F 0x804 -| samtools sort - >{bam}
        samtools index {bam}
        samtools depth -aa {bam}>{depth_file}
        echo alignment done.
        """
    else:
        alignDB_order = f"""
        outdir={args["o"]}/{args["n"]}
        sample={args["n"]}
        ref={ref}
        minimap2 -t {args["j"]} {minimap_para} -E 8,4 -p 0.1 -N 100000 -a $ref {sub_fastq} > {sam}
        # bwa index $ref
        # bwa mem -R '@RG\\tID:foo\\tSM:bar' -a -t {args["j"]} $ref {sub_fastq} > {sam}
        samtools view -bS -F 0x804  {sam} | samtools sort - >{bam}
        samtools index {bam}
        samtools depth -aa {bam}>{depth_file}
        rm {sam}
        echo alignment done.
        """

    # print (alignDB_order)
    # if the depth_file is not detected or empty, then run the alignment
    if not os.path.exists(depth_file) or not os.path.exists(bam) or os.path.getsize(depth_file) == 0 or os.path.getsize(bam) == 0:
        os.system(alignDB_order)
    else:
        print("Depth and BAM file is detected.")
    # os.system(alignDB_order)

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

def read_bin_map2db(args, my_db):
    # map raw fastq to the database of all alleles
    read_type = Read_Type(args["seq_tech"], args["y"], args["RNA_type"])
    minimap_para = read_type.get_minimap2_param()
    minimap_db = my_db.full_db

    # if args["minimap_index"] == 1 and args["seq_tech"] != 'rna':
    #     ref_index = my_db.full_db[:-5] + args["y"] + ".mmi"
    #     # print ("search the reference index:", ref_index)
    #     if not os.path.isfile(ref_index):
    #         print ("start build Minimap2 index for the reference...")
    #         os.system(f"minimap2 {minimap_para} -d {ref_index} {my_db.full_db} ")
    #     else:
    #         print (f"Detect Minimap2 index for the reference: {ref_index}")
    #     minimap_db = ref_index

    outbam = f"""{args["o"]}/{args["n"]}/{args["n"]}.db.bam"""
    # map raw reads to database
    # if args["seq_tech"] != 'rna':
        # alignDB_order = f"""
        # minimap2 -t {args["j"]} {minimap_para} -a {minimap_db} {args["r"]} |samtools view -bS -o {outbam}
        # echo alignment done.
        # """
    alignDB_order = f"""
    bwa mem -t {args["j"]} {my_db.full_db} {args["r"]} |samtools view -bS -o {outbam}
    echo alignment done.
    """
    # else:
    #     cds_bwa_db = my_db.full_cds_db
    #     alignDB_order = f"""
    #     bwa mem -t {args["j"]} {cds_bwa_db} {args["r"]} |samtools view -bS -o {outbam}
    #     echo alignment done.
    #     """
    # print (alignDB_order)
    os.system(alignDB_order)