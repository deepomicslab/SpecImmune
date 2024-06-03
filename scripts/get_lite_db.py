import pysam


def convert_field_for_allele(allele, digit=8):
    array = allele.split(":")
    # print ("allele", allele, array)
    if len(array) < 2:
        return allele
    one = array[0]
    two  = array[0] + ":" + array[1]
    if len(array) >= 3:
        three  = array[0] + ":" + array[1]+ ":" + array[2]
    else:
        three  = two

    if digit == 2:
        allele = one
    elif digit == 4:
        allele = two
    elif digit == 6:
        allele = three
    # mylist[j] = three
    elif digit == 8:
        allele = allele   
    return allele

def extract_contigs(input_fasta, output_fasta):
    """
    Extract contigs with specified names from a FASTA file and save them to a new file.

    Args:
        input_fasta (str): Path to the input FASTA file.
        output_fasta (str): Path to the output FASTA file.
        contig_names (list): List of contig names to extract.

    Returns:
        None
    """
    contig_names = {}
    with pysam.FastaFile(input_fasta) as fasta_in, open(output_fasta, 'w') as fasta_out:
        for contig in fasta_in.references:
            lite_contig = convert_field_for_allele(contig, 4)
            if contig in contig_names:
                continue
            seq = fasta_in.fetch(contig)
            fasta_out.write(f'>{contig}\n{seq}\n')
            contig_names[lite_contig] = 1


            
if __name__ == "__main__":
    # Example usage
    input_fasta = '/mnt/d/HLAPro_backup/Nanopore_optimize/SpecComplex/db/HLA/ref/HLA.extend.fasta'
    output_fasta = '/mnt/d/HLAPro_backup/Nanopore_optimize/SpecComplex/db/HLA/ref/HLA.select.lite.fasta'


    extract_contigs(input_fasta, output_fasta)