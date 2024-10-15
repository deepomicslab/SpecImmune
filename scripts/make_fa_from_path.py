import argparse, re
from collections import OrderedDict
import pysam


def parse_segs():
    segs_dict = OrderedDict()
    with open(args.segs, 'r') as segs_f:
        for line in segs_f:
            line = line.strip().split('\t')
            chrom, pos1, pos2, svtype, svlen = line
            seg_str = "{}-{}".format(pos1, pos2)
            segs_dict[seg_str] = svtype
    return segs_dict, chrom

def reverse_complement(sequence):
    sequence = sequence.upper()  # Convert input sequence to uppercase
    complement = {
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C'
    }
    return ''.join([complement[base] for base in sequence[::-1]])

def parse_path():
    path_order=[]
    with open(args.path, 'r') as path_f:
        lines = path_f.readlines()
    lines.reverse()
    for idx, line in enumerate(lines):
        print(line.strip().split())
        if idx == 0:
            for seg in line.strip().split():
                path_order.append(seg)
    return path_order

def make_fa(path_order):
    out_fa = open(args.out_fa, 'w')
    seg_fa = pysam.FastaFile(args.contigs)
    out_fa.write('>{}\n'.format("{}".format(chrom)))
    all_fasta = ""
    for seg in path_order:
        seg_name, seg_ori = seg[:-1], seg[-1]
        contig_name = f"{chrom}:"+seg_name
        contig_fasta = seg_fa.fetch(contig_name) if seg_ori == '+' else reverse_complement(seg_fa.fetch(contig_name))
        all_fasta += contig_fasta
    print("first 1000 bp: {}".format(all_fasta[:1000]))
    out_fa.write(all_fasta)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Filter BND events not in the same chromosome in a SV VCF file')
    parser.add_argument('path', help='Input VCF file')
    parser.add_argument('segs', help='Input VCF file')
    parser.add_argument('contigs', help='Input VCF file')
    parser.add_argument('out_fa', help='Input VCF file')
    parser.add_argument('hap_idx', help='Input VCF file')

    args = parser.parse_args()

    parse_segs()
    seg_dict, chrom = parse_segs()
    path_order = parse_path()
    make_fa(path_order)