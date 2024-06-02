import os
## canu -nanopore  A.long_read.fq.gz -d test -p test genomeSize=4000



def assembly(fq, assmbly_dir, prefix, genome_size, seq_type="nanopore"):

    """
    -pacbio      <files>
    -nanopore    <files>
    -pacbio-hifi <files>
    """

    cmd = f"""
        canu -{seq_type}  {fq} -d {assmbly_dir} -p {prefix} genomeSize={genome_size}

        echo result is in {assmbly_dir}/{prefix}.contigs.fasta
    """
    contig = f"{assmbly_dir}/{prefix}.contigs.fasta"
    

    os.system(cmd)

    return contig


# for each fastq in a folder, run the assembly by the function:assembly with the following parameters
def assembly_all(fq_dir, assmbly_dir, genome_size, seq_type="nanopore"):
    for fq in os.listdir(fq_dir):
        if fq.endswith(".fq.gz") or fq.endswith(".fastq.gz"):
            prefix = fq.split(".")[0]
            fq = os.path.join(fq_dir, fq)
            assembly(fq, assmbly_dir, prefix, genome_size, seq_type)