import pysam
import os
import argparse
import re 
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

parser=argparse.ArgumentParser()
parser.add_argument("-i", required=True, help="input bamfile")
parser.add_argument("-o", required=True, help="Output bamfile")

def re_mapping_qual_reads(inbam, outbam):
    nf = pysam.AlignmentFile(inbam, 'rb')
    outf = pysam.AlignmentFile(outbam, "wb", header=nf.header) 
    
    for r in nf:
        if r.is_supplementary:
            continue
        if r.is_unmapped:
            continue
        if r.is_read1:
            tag=1
        if r.is_read2:
            tag=2
        if r.mapping_quality < 30:
            r.mapping_quality = 30

        if r.reference_name != r.next_reference_name:
            continue

        outf.write(r)

    nf.close()
    outf.close()



if __name__ == "__main__":
    args=parser.parse_args()
    if args.i:
        inbam = args.i
    if args.o:
        outbam = args.o
    re_mapping_qual_reads(inbam, outbam)


