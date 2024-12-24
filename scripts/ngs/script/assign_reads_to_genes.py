"""
assin reads to original gene loci with the alignment 
"""

import sys
import os
import pysam
import numpy as np
import re
import time
import gzip

B = ['A', 'C', 'G', 'T']
BN = ['A', 'C', 'G', 'T', 'N']

def get_names():
    with open(reads_file, 'r') as infile:
        n = infile.read().splitlines()
    if '' in n:
        n.remove('')
    return n

def count_alignment(alignment):
    match_num = 0
    soft_num = 0
    all_num = 0
    for ci in alignment.cigar:
        if ci[0] == 0:
            match_num += ci[1]
        elif ci[0] == 4:
            soft_num += ci[1]
        all_num += ci[1]

    mis_NM = 0
    for ta in alignment.get_tags():
        if ta[0] == 'NM':
            match_num -= ta[1]
            mis_NM += ta[1]
    focus_len = all_num - soft_num

    return mis_NM, soft_num, match_num, focus_len

def check_score(dict, options, name, pair_dict):
    gene_dict = {}
    for align in dict.keys():
        # print (l)
        if pair_dict[align] < 2:
            dict[align] = 0
        gene = align.split('*')[0]
        if gene == "KIR2DL5A" or gene == "KIR2DL5B":
            gene = "KIR2DL5"
        if gene not in gene_dict.keys():
            gene_dict[gene] = dict[align]
        else:
            if gene_dict[gene] < dict[align]:
                gene_dict[gene] = dict[align]
    new_l = sorted(gene_dict.items(), key=lambda gene_dict:gene_dict[1], reverse = True)
    if float(new_l[0][1]) < 0.1:
        print (name, 'mismatch', new_l)
        return 'REMOVE'
    elif len(new_l) == 1:
        print (name, 'other', new_l)
        return new_l[0][0]
    elif new_l[0][1] - new_l[1][1] < options.diff_score:
        print (name, 'diff', new_l)
        return 'REMOVE'
    else:
        print (name, 'other', new_l)
        return new_l[0][0]

class Each_read():
    """
    for each pair-end reads
    record all the alignment scores
    assign the pair-end read to gene loci by alignment scores
    """

    def __init__(self):
        self.dict = {}
        self.pair_dict = {}
        self.len_dict = {}    
        self.read_name = ''   

    def add_one_alignment(self, alignment):
        """
        given a alignment record,
        save the mapping situations to the read's dicts.
        """
        flag = True
        if alignment.is_unmapped or alignment.reference_name != alignment.next_reference_name:
            flag = False
        t_name = alignment.reference_name
        gene = "NA"
        if  t_name == None:
            gene = "NA"
        else:
            gene = t_name.split('*')[0]
        #print(t_name, alignment.query_name, gene)
        mis_NM, soft_num, match_num, focus_len = count_alignment(alignment)
        print(t_name, alignment.query_name,mis_NM,soft_num,match_num,focus_len)
        if soft_num > 0 and options.soft > 0:
            flag = False
        if mis_NM > options.max_nm:
            flag = False
        if mis_NM > 1 and gene == "CYP2D6":
            flag = False
        #if mis_NM < 50 and gene == "CYP8A1":
        #    flag = True
        #if mis_NM < 50 and gene == "CYP3A5":
        #    flag = True
        #if mis_NM < 50 and gene == "CYP4B1":
        #    flag = True
        if mis_NM > 3 and gene == "CYP4F2":
            flag = False
        if mis_NM > 1 and gene == "CYP2C9":
            flag = False
        if mis_NM > 1 and gene == "CYP2C19":
            flag = False
        if mis_NM > 1 and gene == "CYP2A6":
            flag = False
        if mis_NM > 2 and gene == "CYP2B6":
            flag = False
        if flag:
            self.read_name = alignment.query_name
            if t_name not in self.dict.keys():
                self.dict[t_name] = match_num#round(s,3)
                self.len_dict[t_name] = focus_len
                self.pair_dict[t_name] = 1
            else:
                self.dict[t_name] += match_num#round(s,3)
                self.len_dict[t_name] += focus_len
                self.pair_dict[t_name] += 1
    
    def assign(self):
        """
        determine which locus the read should be assigned
        """
        
        flag = True
        if len(self.dict) == 0:
            flag = False
        for key in self.dict.keys():
            # print (len_dict[key])
            if self.len_dict[key] <= 0:  #make sure the reads is paired mapped.
                self.dict[key] = 0
            else:
                self.dict[key] = float(self.dict[key])/self.len_dict[key]
        if flag:
            first_align = check_score(self.dict, options, self.read_name, self.pair_dict)
        else:
            first_align = 'REMOVE'
        return first_align

def assign_fastq(file, gene, index, assign_dict):
    """
    extract gene-specific reads from the raw read file,
    generate gene-specific fastq files
    read name in the original fastq should be
    @<read name> or @<read name>/1
    """
    i = 0
    #gene = 'A'
    outfile = options.outdir + '/%s.R%s.fq'%(gene, index)
    out = open(outfile, 'w')
    flag = False
    if file.split(".")[-1] == "gz":
        f = gzip.open(file,'rt')
    else:
        f = open(file)
    for line in f:
        line = line.strip()
        if i % 4 == 0:
            if re.search('/1',line) or re.search('/2',line):
                read_name = line.split()[0][1:-2]
            else:
                read_name = line.split()[0][1:]
            if read_name in assign_dict.keys() and assign_dict[read_name] == gene:
                flag = True
                num = 1
                print (line, file = out)
        elif flag:
            print (line, file = out)
            num += 1
            if num == 4:
                flag = False
        i += 1
    out.close()
    os.system('gzip -f %s'%(outfile))

def main():
    print ('start assigning reads...')
    read_dict = {} # record the aligment scores for each read
    assign_dict = {} # record assigned gene for each read
    bamfile = pysam.AlignmentFile(options.bam, 'rb')
    t0 = time.time()

    # record alignment scores for all the read
    for alignment in bamfile:
        t_name = alignment.reference_name
        alignment.query_name = alignment.query_name.split("/")[0]
        read_name = alignment.query_name
        # print (read_name)
        if read_name not in read_dict.keys():
            
            new_read = Each_read()
            read_dict[read_name] = new_read
        read_dict[read_name].add_one_alignment(alignment)
    bamfile.close()
    t1 = time.time()
    print ("read bam cost %s"%(t1 - t0))

    # assign genes for each read
    for read_name in read_dict:
        assigned_locus = read_dict[read_name].assign()
        if assigned_locus == 'REMOVE':
            continue
        assign_dict[read_name] = assigned_locus
    if options.gene_class == "HLA":
        Genes_list = [ 'HLA-A', 'HLA-B', 'HLA-C', 'HLA-DMA', 'HLA-DMB', 'HLA-DOA', 'HLA-DOB', 'HLA-DPA1', 'HLA-DPB1', 'HLA-DPB2', 'HLA-DQA1', 'HLA-DQB1', 'HLA-DRA', 'HLA-DRB1', 'HLA-DRB3', 'HLA-DRB4', 'HLA-DRB5', 'HLA-E', 'HLA-F', 'HLA-G', 'HLA-H', 'HLA-J', 'HLA-K', 'HLA-L', 'HLA-P', 'HLA-V', 'HLA-DQA2', 'HLA-DPA2', 'HLA-N', 'HLA-S', 'HLA-T', 'HLA-U', 'HLA-W', 'MICA', 'MICB', 'TAP1', 'TAP2', 'HFE' ]
    if options.gene_class == "CYP":
        Genes_list = [ 'CYP19A1', 'CYP1A1', 'CYP1B1', 'CYP26A1', 'CYP2A13', 'CYP2A6', 'CYP2B6', 'CYP2C19', 'CYP2C8', 'CYP2C9', 'CYP2D6', 'CYP2F1', 'CYP2J2', 'CYP2R1', 'CYP2S1', 'CYP2W1', 'CYP3A4', 'CYP3A43', 'CYP4A22', 'CYP4B1', 'CYP4F2', 'CYP8A1', 'CYP3A5', 'CYP3A7' ]
    if options.gene_class == "KIR":
        Genes_list = ['KIR2DL1', 'KIR2DL2', 'KIR2DL3', 'KIR2DL4', 'KIR2DL5', 'KIR2DP1', 'KIR2DS1', 'KIR2DS2', 'KIR2DS3', 'KIR2DS4', 'KIR2DS5', 'KIR3DL1', 'KIR3DL2', 'KIR3DL3', 'KIR3DP1', 'KIR3DS1']

    # generate gene-specific fastq
    for gene in Genes_list:
        assign_fastq(options.fq1, gene, 1, assign_dict)
        assign_fastq(options.fq2, gene, 2, assign_dict)
    t2 = time.time()
    print ("read assigment cost %s"%(t2 - t0))

if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(description='Read assignment')
    parser.add_argument('-a', '--gene_class', help='gene class', required=True)
    parser.add_argument('-b', '--bam', help='bam file', required=True)
    parser.add_argument('-o', '--outdir', help='outdir', required=True)  
    parser.add_argument('-1', '--fq1', help='bin dir', required=True) 
    parser.add_argument('-2', '--fq2', help='bin dir', required=True)  
    parser.add_argument('-n', '--bin_dir', help='bin dir', required=True)  
    parser.add_argument('-s', '--soft', help='remove soft-clip reads', required=False, default =1,  type=int)
    parser.add_argument('-nm', '--max_nm', help='MAX NM', required=False, default = 2, type=int)
    parser.add_argument('-d', '--diff_score', help='The score for the best gene must be at least this higher than the second gene', required=False, default = 0.5, type=float)
    options = parser.parse_args()

    main()
