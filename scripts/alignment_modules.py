
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


if __name__ == "__main__":

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