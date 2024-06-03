


class Alignment():

    def __init__(self, args):
                self.args = args

        self.reference = 

        self.bamfile = 
        self.sort_bam = None
        self.threads = args["j"]
        self.tool = tool



    def index_db(self):
        ref_index = self.reference[:-5] + self.args["y"] + ".mmi"
        # print ("search the reference index:", ref_index)
        if not os.path.isfile(ref_index):
            print ("start build Minimap2 index for the reference...")
            os.system(f"minimap2 {minimap_para} -d {ref_index} {self.db} ")
        else:
            print (f"Detect Minimap2 index for the reference: {ref_index}")
        self.db = ref_index

    

    def run(self):
    
        if args["minimap_index"] == 1:
            self.index_db()
        # map raw reads to database
        alignDB_order = f"""
        fq={parameter.raw_fq}
        ref={self.db}
        outdir={parameter.outdir}
        sample={parameter.sample}
        # minimap2 -t {self.threads} {minimap_para} -a $ref $fq |samtools view -bS -o {self.bamfile}
        bwa mem -R '@RG\\tID:foo\\tSM:bar' -t {self.threads} {my_db.full_db} $fq |samtools view -bS -o {self.bamfile}
        echo alignment done.
        """
        # print (alignDB_order)
        os.system(alignDB_order)

    def sort_bam(self):

    
    def count_depth(self):

        cmd = """
        sample=%s
        outdir=%s
        hla=%s
        hla_ref=%s
        minimap2 -t %s %s -a $hla_ref $outdir/$hla.%s.fq.gz | samtools view -bS -F 0x800 -| samtools sort - >$outdir/$hla.bam
        samtools index $outdir/$hla.bam
        samtools depth -d 1000000 -aa $outdir/$hla.bam >$outdir/$hla.depth
        """%(parameter.sample, parameter.outdir, gene, my_db.get_gene_alleles(gene), parameter.threads, minimap_para, args["a"])
        os.system(cmd)


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


reference = "reference.fasta"  # Replace with the path to your reference FASTA file
fastq = "reads.fastq"  # Replace with the path to your raw FASTQ file

aligner = LongReadAligner(reference, fastq, alignment_method="minimap2")

# Index the reference (if not already indexed)
aligner.index_reference()

# Perform alignment
output_sam = "aligned_reads.sam"
aligner.align_reads(output_sam)

# Sort the resulting BAM file
input_bam = "aligned_reads.sam"
output_bam = "sorted_reads.bam"
aligner.sort_bam(input_bam, output_bam)

# Calculate depth from the sorted BAM file
output_depth = "depth.txt"
aligner.calculate_depth(output_bam, output_depth)

# Print the depth file
with open(output_depth, 'r') as depth_file:
    for line in depth_file:
        print(line.strip())