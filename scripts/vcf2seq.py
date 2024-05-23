import argparse
import pysam
from Bio.Seq import Seq

def reverse_complement(seq):
    """ Return the reverse complement of the input DNA sequence. """
    return str(Seq(seq).reverse_complement())

def process_vcf(input_vcf, ref_fasta, output_vcf):
    vcf_in = pysam.VariantFile(input_vcf)
    fasta = pysam.FastaFile(ref_fasta)
    vcf_out = pysam.VariantFile(output_vcf, 'w', header=vcf_in.header)

    for record in vcf_in.fetch():
        # Replace 'N' with the actual base from the reference genome
        if record.ref == 'N':
            actual_ref = fasta.fetch(record.contig, record.start, record.start + 1)
            record.ref = actual_ref
        else:
            actual_ref = record.ref

        # Process alt alleles
        new_alts = []
        for alt in record.alts:
            if alt == "<INV>":
                # continue
                # Get the reverse complement of the actual reference sequence
                actual_ref_seq = fasta.fetch(record.contig, record.start, record.stop)
                #record.ref = actual_ref_seq
                new_alts.append(reverse_complement(actual_ref_seq))
            elif alt in ("<INS>", "<DUP>", "<BND>"):
                # Skip these variant types
                continue
            else:
                new_alts.append(alt)

        # Update alts if there are valid alternatives left
        if new_alts:
            record.alts = tuple(new_alts)
            vcf_out.write(record)

    vcf_in.close()
    vcf_out.close()
    fasta.close()

def main():
    parser = argparse.ArgumentParser(description="Process a VCF file to adjust REF bases and handle specific variant types, including correct handling of <INV> variants.")
    parser.add_argument("input_vcf", help="Input VCF file path")
    parser.add_argument("ref_fasta", help="Reference genome FASTA file path")
    parser.add_argument("output_vcf", help="Output VCF file path")

    args = parser.parse_args()

    process_vcf(args.input_vcf, args.ref_fasta, args.output_vcf)

if __name__ == "__main__":
    main()
