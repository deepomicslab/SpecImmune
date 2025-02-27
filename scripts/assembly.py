import os
from Bio import SeqIO
import argparse
import sys
import pysam

from alignment_modules import Read_Type
from db_objects import My_db

## canu -nanopore  A.long_read.fq.gz -d test -p test genomeSize=4000

## given a fasta file, only retain the first contig using biopython

def get_first_contig(fasta_file, output_file):
    with open(output_file, "w") as out:
        for record in SeqIO.parse(fasta_file, "fasta"):
            out.write(f">{record.id}\n{record.seq}\n")
            break
        out.close()

def get_longest_read_length(fq):
    longest_read = 0
    with pysam.FastxFile(fq) as fh:
        for entry in fh:
            if len(entry.sequence) > longest_read:
                longest_read = len(entry.sequence)
    return longest_read


def assembly(fq, outdir, genome_size, gene, index, flye_param, allele):
    prefix = f"{gene}_{index}"
    assmbly_dir = outdir + f"/{prefix}"
    # cmd = f"""
    #     canu -{seq_type}  {fq} -d {assmbly_dir} -p {prefix} genomeSize={genome_size}

    # """
    # contig = f"{assmbly_dir}/{prefix}.contigs.fasta"
    longest_read = get_longest_read_length(fq)

    if longest_read > 1000:

        cmd = f"""
            flye --min-overlap 1000 {flye_param} {fq} --out-dir {assmbly_dir} --genome-size {genome_size} --threads {args["j"]}
        """

        print (cmd)
        os.system(cmd)
        contig = f"{assmbly_dir}/assembly.fasta"

        ### only keep the first contig in the final result if there are multiple contigs
        get_first_contig(contig, f"{outdir}/hla.allele.{index}.{gene}.fasta")
    else:
        cmd = f"""
            samtools faidx {my_db.full_db}  {allele} > {outdir}/hla.allele.{index}.{gene}.fasta
        """
        print (cmd)
        os.system(cmd)

    return contig


# for each fastq in a folder, run the assembly by the function:assembly with the following parameters
def assembly_all(outdir, genome_size, flye_param):
    gene_index_dict = {}
    for fq in os.listdir(outdir):
        if fq.endswith(".fq.gz"):
            allele = fq.split(".")[0]
            if len(allele.split("*")) > 1:
                gene = allele.split("*")[0]
            else:
                continue
            if gene in gene_index_dict:
                gene_index_dict[gene] = 2
            else:
                gene_index_dict[gene] = 1
            fq = os.path.join(outdir, fq)
            
            assembly(fq, outdir, genome_size, gene, gene_index_dict[gene], flye_param, allele)

def assembly_single(fq, outdir, genome_size, flye_param, asm_dir):
    # cmd = f"""
    #     canu -{seq_type}  {fq} -d {assmbly_dir} -p {prefix} genomeSize={genome_size}

    # """
    # contig = f"{assmbly_dir}/{prefix}.contigs.fasta"
    longest_read = get_longest_read_length(fq)

    if longest_read > 1000:

        cmd = f"""
            flye --min-overlap 1000 {flye_param} {fq} --out-dir {asm_dir} --genome-size {genome_size} --threads {args["j"]}
        """

        print (cmd)
        os.system(cmd)
        contig = f"{asm_dir}/assembly.fasta"

        ### only keep the first contig in the final result if there are multiple contigs
        # get_first_contig(contig, f"{outdir}/hla.allele.{index}.{gene}.fasta")
    

    return contig


def assembly_from_whatshap(outdir, genome_size, flye_param):
    gene_index_dict = {}
    h0_fq=outdir + f"/{args['n']}.dv.h0.fq.gz"
    asm_prefix = outdir + f"/{args['n']}_h0"
    assembly_single(h0_fq, outdir, genome_size, flye_param, asm_prefix)
    h1_fq=outdir + f"/{args['n']}.dv.h1.fq.gz"
    asm_prefix = outdir + f"/{args['n']}_h1"
    assembly_single(h1_fq, outdir, genome_size, flye_param, asm_prefix)



if __name__ == "__main__":  
    parser = argparse.ArgumentParser(description="assembly.", add_help=False, \
    usage="python3 %(prog)s -h", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    required = parser.add_argument_group("Required arguments")
    optional = parser.add_argument_group("Optional arguments")

    required.add_argument("-o", type=str, help="The output folder to store the typing results.", metavar="\b", default="./output")
    required.add_argument("-n", type=str, help="Sample ID", metavar="\b")
    required.add_argument("-i", type=str, help="HLA,KIR,CYP",metavar="\b", default="HLA")
    optional.add_argument("-j", type=int, help="Number of threads.", metavar="\b", default=5)
    optional.add_argument("-y", type=str, help="Read type, [nanopore|pacbio].", metavar="\b", default="pacbio")
    optional.add_argument("--seq_tech", type=str, help="Amplicon sequencing or WGS sequencing [wgs|amplicon].", metavar="\b", default="wgs")
    optional.add_argument("-rt", "--RNA_type", type=str, help="traditional,2D,Direct,SIRV",metavar="\b", default="traditional")
    optional.add_argument("--db", type=str, help="db dir.", metavar="\b", default=sys.path[0] + "/../db/")
    optional.add_argument("--hg38", type=str, help="referece fasta file, used by IG_TR typing, generated by extract_VDJ_segments_from_hg38.py", metavar="\b", \
                          default=sys.path[0] + "/../VDJ_ref/IG_TR.segment.fa")
    optional.add_argument("-h", "--help", action="help")
    args = vars(parser.parse_args()) 

    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(0)

    read_type = Read_Type(args["seq_tech"], args["y"],  args["RNA_type"])
    my_db = My_db(args)
    flye_param = read_type.get_flye_param()
    outdir = args["o"] + "/" + args["n"]

    # assembly_all("/mnt/d/HLAPro_backup/Nanopore_optimize/data/sim_hap/reads/pacbio_dp50_acc95_1", 10000, flye_param)
    # assembly_all(outdir, 10000, flye_param)
    assembly_from_whatshap(outdir, 10000, flye_param)


