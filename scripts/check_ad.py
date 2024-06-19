import pysam
import sys

def calculate_haplotype_frequencies(vcf_file):
    haplotype1_support = 0
    haplotype2_support = 0
    total_heterozygous_sites = 0

    vcf = pysam.VariantFile(vcf_file)

    for record in vcf:
        for sample in record.samples:
            gt = record.samples[sample]['GT']
            ad = record.samples[sample]['AD']
            
            # Check if genotype is heterozygous
            if gt[0]!=gt[1]:
                total_heterozygous_sites += 1
                
                # Compare allele depths to determine haplotype support
                if ad[gt[0]] > ad[gt[1]]:
                    haplotype1_support += 1
                    #print(record.pos, gt, ad, "sup 0")
                else:
                    haplotype2_support += 1
                    #print(record.pos, gt, ad, "sup 1")

    if total_heterozygous_sites == 0:
        print("No heterozygous sites found.")
        return

    haplotype1_frequency = haplotype1_support / total_heterozygous_sites
    haplotype2_frequency = haplotype2_support / total_heterozygous_sites
    # Determine the most supported allele index
    if haplotype1_support > haplotype2_support:
        most_supported_allele_index = 1
    else:
        most_supported_allele_index = 2

    print(most_supported_allele_index)
    #print(f"Haplotype1 Frequency: {haplotype1_frequency:.2f}")
    #print(f"Haplotype2 Frequency: {haplotype2_frequency:.2f}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <vcf_file>")
        sys.exit(1)
    
    vcf_file = sys.argv[1]
    calculate_haplotype_frequencies(vcf_file)