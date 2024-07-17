###

"""
define a func to load this file:
#version=pharmvar-6.1.2.1
Haplotype Name  Gene    rsID    ReferenceSequence       Variant Start   Variant Stop    Reference Allele        Variant Allele  Type
CYP2D6*1        CYP2D6          REFERENCE       .
CYP2D6*1.001    CYP2D6          REFERENCE       .
CYP2D6*1.002    CYP2D6  rs28371732      NC_000022.11    42126963        42126963        C       T       substitution
CYP2D6*1.003    CYP2D6  rs150163869     NC_000022.11    42128813        42128813        G       A       substitution
CYP2D6*1.004    CYP2D6  rs28371718      NC_000022.11    42128216        42128216        G       T       substitution
CYP2D6*1.005    CYP2D6  rs111606937     NC_000022.11    42128922        42128922        A       G       substitution
CYP2D6*1.006    CYP2D6  rs28371699      NC_000022.11    42130482        42130482        C       A       substitution
CYP2D6*1.006    CYP2D6  rs28371702      NC_000022.11    42129950        42129950        A       C       substitution
"""
def load_pharmvar(file):
    with open(file, 'r') as f:
        f.readline()
        for line in f:
            line = line.strip().split()
            print (line)
            haplotype_name = line[0]
            gene = line[1]
            rsID = line[2]
            if rsID == "REFERENCE":
                pass
            else:
                ref_seq = line[3]
                variant_start = line[4]
                variant_stop = line[5]
                ref_allele = line[6]
                variant_allele = line[7]
                variant_type = line[8]
            print(haplotype_name, gene, rsID, ref_seq, variant_start, variant_stop, ref_allele, variant_allele, variant_type)

file = "/mnt/d/HLAPro_backup/Nanopore_optimize/data/pharmvar-6.1.2.1/CYP2D6/GRCh38/CYP2D6.NC_000022.11.haplotypes.tsv"
load_pharmvar(file)