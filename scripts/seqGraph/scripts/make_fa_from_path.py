import pysam
import sys


def rev_comp(s):
    sc = ''
    for i in reversed(s):
        if i == 'A':
            sc += 'T'
        elif i == 'C':
            sc += 'G'
        elif i == 'G':
            sc += 'C'
        elif i == 'T':
            sc += 'A'
    return sc

def get_len_long_substr(str1, str2):
    substring = ''
    len_str1 = len(str1)
    if len_str1 > 0:
        for i in range(len_str1):
            for j in range(len_str1 - i + 1):
                if j > len(substring) and all(str1[i:i + j] in x for x in [str1, str2]):
                    substring = str1[i:i + j]

        idx1 = str1.index(substring)
        idx2 = str2.index(substring)

    return idx1, idx2, len(substring), substring

def negative_n(contig1,contig2, fa,ncount):
    c1 = fa.fetch(contig1[:-1])
    if contig1[-1] == "-":
        c1 = rev_comp(fa.fetch(contig2[:-1]))
    c2 = fa.fetch(contig2[:-1])
    if contig2[-1] == "-":
        c2 = rev_comp(fa.fetch(contig2[:-1]))

    c1_sub = c1[len(c1) - (ncount+5):]
    c2_sub = c2[0:ncount+5]

    c1_idx,c2_idx,len_sub,_2= get_len_long_substr(c1_sub, c2_sub)
    if abs(len_sub, ncount) > 5:
        return -1
    else:
        return ncount - c1_idx


def parser_gap(gapf):
    r = {}
    gapin = open(gapf)
    for line in gapin.readlines():
        arr = line.strip().split("\t")
        r[arr[0]+arr[1]] = int(arr[2])
    return r
if __name__ == "__main__":
    fa = pysam.FastaFile(sys.argv[1])
    res = open(sys.argv[4], "w")
    paths = open(sys.argv[2], 'r')
    gaps = sys.argv[3]
    # fa = pysam.FastaFile("/Users/troye/Documents/Recycler/bin/sim800/all.fa")
    # res = open("/Users/troye/Documents/Recycler/bin/sim800/neo800.fa", "w")
    # paths = open("/Users/troye/Documents/tmp/seqGraph/r.debug.c.522446.test.txt", 'rb')
    gaps_r = parser_gap(gaps)
    Ns_count = {}

    for i, line in enumerate(paths.readlines()):
        if line.startswith("iter") or line.startswith("self"):
            continue
        contigs = line.strip().split('\t')
        for idx in range(len(contigs) -1 ):
            key = contigs[idx] + contigs[idx + 1]
            if i not in Ns_count.keys():
                Ns_count[i] = []
            if key in gaps_r.keys():
                Ns_count[i].append(gaps_r[key])
            else:
                Ns_count[i].append(100)
        if i not in Ns_count.keys():
            Ns_count[i] = []
        Ns_count[i].append(0)

    paths.close()
    for idx, line in enumerate(open(sys.argv[2]).readlines()):
        if line.startswith("iter") or line.startswith("self"):
            continue
        fasta = ""
        contigs = line.strip().split('\t')
        clen = len(contigs)
        ns = Ns_count[idx]
        for cidx,contig in enumerate(contigs):
            # print(contig)
            if contig[-1] == "+":
                try:
                    fasta += fa.fetch(contig[:-1])
                    if cidx != clen - 1:
                        if (ns[0] < 0):
                            next_contig = contigs[cidx + 1]
                            cut_len = negative_n(contig,next_contig, fa, ns[0])
                            if cut_len > 0:
                                fasta = fasta [0:len(fasta) - cut_len]
                        else:
                            fasta += 'N'*ns[0]
                except Exception as e:
                    name = "_".join(contig[:-1].split('_')[:-1])
                    print("not found: ", name)
                    fasta += fa.fetch(name)
                    if cidx != clen - 1:
                        fasta += 'N'*ns[0]
            else:
                try:
                    fasta += rev_comp(fa.fetch(contig[:-1]))
                    if cidx != clen - 1:
                        fasta += 'N'*ns[0]
                except Exception as e:
                    name = "_".join(contig[:-1].split('_')[:-1])
                    print("not found: ", name)
                    fasta += rev_comp(fa.fetch(name))
                    if cidx != clen - 1:
                        if (ns[0] < 0):
                            next_contig = contigs[cidx + 1]
                            cut_len = negative_n(contig,next_contig, fa, ns[0])
                            if cut_len > 0:
                                fasta = fasta [0:len(fasta) - cut_len]
                        else:
                            fasta += 'N'*ns[0]
            ns.pop(0)
        res.write(">res_" + str(idx + 1) + "_" + str(len(fasta)) + "\n")
        res.write(fasta + "\n")

    #
    # i = 0
    # for line in paths.readlines():
    #     fasta = ""
    #     if i % 2 != 0:
    #         contigs = line.strip().split('\t')
    #         for contig in contigs:
    #             # print(contig)
    #             if contig[-1]=="+":
    #                 print(i, contig[:-1], contig[-1])
    #                 try:
    #                     fasta += fa.fetch(contig[:-1])
    #                 except Exception as e:
    #                     name = "_".join(contig[:-1].split('_')[:-1])
    #                     print(i, name)
    #                     fasta += fa.fetch(name)
    #             else:
    #                 print(i, contig[:-1], contig[-1])
    #                 try:
    #                     fasta += rev_comp(fa.fetch(contig[:-1]))
    #                 except Exception as e:
    #                     name = "_".join(contig[:-1].split('_')[:-1])
    #                     print(i, name)
    #                     fasta += rev_comp(fa.fetch(name))
    #         res.write(">res"+str(i+1)+"_"+str(len(fasta))+"\n")
    #         res.write(fasta+"\n")
    #     i += 1

    # chx=fa.fetch("EDGE_1298087_length_128_cov_268.000000")
    # print(chx)
    # with pysam.FastxFile("/Users/troye/Documents/Recycler/bin/522446/all.fa") as fin:
    #     chrx= fin.fetch("EDGE_61_length_838_cov_170.571027")
    #     print(chrx)
    # for r in fin:
    #     print(r)
    res.close()
    paths.close()
