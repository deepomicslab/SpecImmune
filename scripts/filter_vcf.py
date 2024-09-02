"""
Given a vcf file, remove the site with  minor allele frequency less than 0.3
"""

from pysam import VariantFile
import numpy as np
import sys

def filter_vcf(vcf_file, output_file, maf=0.3):
    maf = float(maf)
    ## use pysam
    bcf_in = VariantFile(vcf_file)  # auto-detect input format
    bcf_out = VariantFile(output_file, 'w', header=bcf_in.header)
    sample=list(bcf_in.header.samples)[0]
    for rec in bcf_in.fetch():

        geno=rec.samples[sample]['GT']  
        depth=rec.samples[sample]['AD']   
        depth_sum = sum(depth)
        freq = np.array(depth)/depth_sum

        if len(depth) > 2:  ## retail sites with only two 2 alleles
           bcf_out.write(rec)
        elif min(freq) < maf:
            # if 'UG' in rec.samples[sample]:
            genotype = rec.samples[sample]['GT']
            genotype = sorted(genotype)
            # print (freq, genotype)
            if freq[0] >= freq[1]:
                rec.samples[sample]['GT'] = (genotype[0], genotype[0])
            else:
                rec.samples[sample]['GT'] = (genotype[1], genotype[1])
            rec.samples[sample].phased=True
            bcf_out.write(rec)
        else:
            bcf_out.write(rec)

        
        # print (depth, freq)


if __name__ == "__main__":
    # vcf_file = "/mnt/d/HLAPro_backup/Nanopore_optimize/vdj_results_tcr2/NA10831/genes/TRBV6-6.phase.vcf.gz"
    # output_file = "/mnt/d/HLAPro_backup/Nanopore_optimize/vdj_results_tcr2/NA10831/genes/TRBV6-6.phase.filtered.vcf.gz"
    # maf=0.3

    vcf_file = sys.argv[1]
    output_file = sys.argv[2]
    maf = sys.argv[3]

    
    filter_vcf(vcf_file, output_file, maf)


