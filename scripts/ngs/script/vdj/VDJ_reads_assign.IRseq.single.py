import sys
import os
import pysam
import numpy as np
import re
from collections import defaultdict

B = ['A', 'C', 'G', 'T']
BN = ['A', 'C', 'G', 'T', 'N']


def extract_reads_name(options):
    reads_file = '%s/read_names.txt'%(options.outdir)
    order = 'samtools view %s | cut -f1 | sort | uniq > %s'%(options.bam, reads_file)
    os.system(order)

def get_names(names):
    with open(names, 'r') as infile:
        n = infile.read().splitlines()
    if '' in n:
        n.remove('')
    return n

def extract_reads(options, read12, VJgene):
    print ('start assigning reads...')
    reads_file = '%s/read_names.txt'%(options.outdir)
    assign_file = '%s/%s.%s.assign_file.txt'%(options.outdir,read12, VJgene)
    out = open(assign_file, 'w')
    n = get_names(reads_file)
    bamfile = pysam.AlignmentFile(options.bam, 'rb')
    name_indexed = pysam.IndexedReads(bamfile)
    name_indexed.build()
    error, total, remove = 0, 0, 0
    error_set = []
    flist = {}
    readstats = defaultdict(dict)
    for name in n:
        try:
            iterator = name_indexed.find(name)
            ddict = {}
            pair_dict = {}
            len_dict = {}
            name12 = name
            for x in iterator:   
                if x.is_unmapped:
                    continue
                #if x.is_read2 and read12 == "read1":
                #    continue
                #if x.is_read1 and read12 == "read2":
                #    continue
                #if x.is_read1:
                #    name12 = name + "." + "1"
                #else:
                #    name12 = name + "." + "2"
                s = 0
                t_name = x.reference_name
                if t_name.find('V') == -1 and VJgene == "V":
                    continue
                if t_name.find('J') == -1 and VJgene == "J":
                    continue

                match_num = 0
                soft_num = 0
                all_num = 0
                for ci in x.cigar:
                    if ci[0] == 0:
                        match_num += ci[1]
                    elif ci[0] == 4:
                        soft_num += ci[1]
                    all_num += ci[1]

                mis_NM = 0
                pattern = re.compile(r'\d+S\S+S$', re.I)
                scigar = re.findall(pattern, x.cigarstring)
                for ta in x.get_tags():
                    if ta[0] == 'NM':
                        match_num -= ta[1] * options.penalty
                        mis_NM += ta[1]
                if t_name.find('V') > -1 and mis_NM > options.max_nm:
                    continue
                #if t_name.find('V') > -1 and len(scigar) == 1:
                #    continue
                if t_name.find('V') == -1 and mis_NM > 1:
                    continue
                focus_len = all_num - soft_num
                rate = mis_NM / (match_num + mis_NM * options.penalty )
                if rate > options.rate:
                    continue

                if t_name not in ddict.keys():
                    ddict[t_name] = match_num
                    len_dict[t_name] = focus_len
                    pair_dict[t_name] = 1
                else:
                    
                    ddict[t_name] += match_num
                    len_dict[t_name] += focus_len
                    pair_dict[t_name] += 1

                readstats[name12][t_name] = ddict[t_name]

            total += 1
            if len(ddict) == 0:
                continue

            first_align = check_score(ddict, options, name, pair_dict)
            #print(name, first_align, ddict)
            if first_align == 'REMOVE':
                remove += 1
                continue
            print (name, first_align, file = out)
            flist[name] = first_align

        except KeyError:
            pass
    out.close()
    return readstats

def select_VDJ_allele(options, readstats, VJgene):
    flist = {}
    clist = defaultdict(int)
    gene_reads = defaultdict(list)
    assign_file1 = '%s/read.%s.assign_file.txt'%(options.outdir, VJgene)
    #assign_file2 = '%s/read2.%s.assign_file.txt'%(options.outdir, VJgene)
    with open(assign_file1) as f:
        for line in f:
            read, gene = line.strip().split(' ')
            #read12 = read + "." + "1"
            flist[read] = gene
            clist[gene] += 1
            gene_reads[gene].append(read)
    f.close()
   # with open(assign_file2) as f:
    #    for line in f:
    #        read, gene = line.strip().split(' ')
    #        read12 = read + "." + "2"
    #        flist[read12] = gene
    #        clist[gene] += 1
    #        gene_reads[gene].append(read12)
    #f.close()
    vdj_file = '%s/read.%s.file.txt'%(options.outdir, VJgene)
    out = open(vdj_file, 'w')
    rdict = defaultdict(list)
    gdict = defaultdict(int)
    gidict = defaultdict(list)
    for rname in flist:
        gid = flist[rname]
        nids = readstats[rname]
        #if rname.endswith("1"):
        #    nids = readstats1[rname]
        #else:
        #    nids = readstats2[rname]
        maxv = max(nids.values())
        if len(nids) == 1:
            gene = list(nids.keys())[0]
            rdict[rname].append(gene)
            gdict[gene] += 1
            rgene = gene.split('*')[0];
            gidict[rgene].append(gene)
        else:
            for k, v in nids.items():
                if v == maxv:
                    rgene = k.split('*')[0]
                    gidict[rgene].append(k)
                    gdict[k] += 1
                    rdict[rname].append(k)
    
    for vv in gdict.keys():
        print(vv,gdict[vv]) # 'IGKV2D-26*01': 33

    selectvdj = defaultdict(list)
    for vdj in gidict.keys():
        vdjs = list(set(gidict[vdj]))
        total = clist[vdj]
        print(vdj, total)
        greads = gene_reads[vdj]
        tdict = {}
        for i in vdjs:
            tdict[i] = gdict[i]
        gene1 = top_n_scores(1,tdict)[0]
        selectvdj[vdj].append(gene1)
        if len(tdict) >= 2:
            if gdict[gene1] >= 0.80 * total:
                continue
            else:
                #gene2 = top_n_scores(2,tdict)[1]
                gene2 = second_gene(rdict, greads, gene1, total, gdict)
                selectvdj[vdj].append(gene2)

    for svdj in selectvdj.keys():
        print(svdj, ",".join(selectvdj[svdj]), file=out)

    out.close()

def second_gene(rdict, greads, first_gene, total, gdict):
    remains = defaultdict(int)
    for read in greads:
        lists = rdict[read]
        if any(first_gene in s for s in lists):
            continue
        else:
            for gene in lists:
                remains[gene] += 1
    gene2 = top_n_scores(1,remains)[0]
    #print(first_gene,gdict[first_gene],gene2,gdict[gene2],total)
    return gene2
                
def top_n_scores(n, sdict):
    lot = [(k,v) for k, v in sdict.items()] #make list of tuple from scores dict
    nl = []
    ml = []
    while len(lot)> 0:
        mit = max(lot, key=lambda x: x[1])
        ml.append(mit[0])
        nl.append(mit)
        lot.remove(nl[-1])
    return ml[0:n]   

def check_score(dict, options, name, pair_dict):
    gene_dict = {}
    for align in dict.keys():
        gene = align.split('*')[0]
        if gene not in gene_dict.keys():
            gene_dict[gene] = dict[align]
        else:
            if gene_dict[gene] < dict[align]:
                gene_dict[gene] = dict[align]
    new_l = sorted(gene_dict.items(), key=lambda gene_dict:gene_dict[1], reverse = True)

    if float(new_l[0][1]) < 0.8:
        print(name, 'REMOVE', new_l)
        return 'REMOVE'
    if len(new_l) == 1:
        print(name, new_l[0][0], new_l)
        return new_l[0][0]
    elif new_l[0][1] - new_l[1][1] < 0.001:# or new_l[0][1] < options.theta_pm:
        print(name, 'REMOVE', new_l)
        return 'REMOVE'
    else:
        print(name, new_l[0][0], new_l)
        return new_l[0][0]

if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(description='Extract reads by read name from bam file')
    parser.add_argument('-b', '--bam', help='bam file', required=True)
    parser.add_argument('-o', '--outdir', help='outdir', required=True)   
    parser.add_argument('-nm', '--max_nm', help='MAX NM', required=False, default = 1, type=int)
    parser.add_argument('-p', '--penalty', help='penalty of mismatch', required=False, default = 2, type=int)
    parser.add_argument('-r', '--rate', help='max mismatch rate', required=False, default=0.05, type=float)

    options = parser.parse_args()
    extract_reads_name(options)
    
    Vreadstats1 = extract_reads(options, "read", "V")
    #Vreadstats2 = extract_reads(options, "read2", "V")
    Jreadstats1 = extract_reads(options, "read", "J")
    #Jreadstats2 = extract_reads(options, "read2", "J")
    #print(readstats)
    select_VDJ_allele(options, Vreadstats1, "V")
    select_VDJ_allele(options, Jreadstats1, "J")
