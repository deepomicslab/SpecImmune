from Bio import SeqIO
import re
import sys
from Bio.Seq import Seq

fain = sys.argv[1]
orderin = open(sys.argv[2])
orderintmp = open(sys.argv[2])
faout = open(sys.argv[3], "w")
blastin = open(sys.argv[4])
blast_ratio = float(sys.argv[5])
gene_hit = sys.argv[6]
plasscore_file = sys.argv[7]
blast_segs = set()
prev_seg = ""
prev_len = 0
prev_ref = ""
res_count = set()
in_faout= []
def contains_gen(line,genes_seg):
    strip_line = line.strip().replace("+","").replace("-","")
    line_arr = strip_line.split("\t")
    for item in line_arr:
        if item in genes_seg:
            return True
    return False
def line_len(line):
    result_len = 0
    vs = re.split(r'[+-]',line)
    for v in vs:
        if v == "":
            continue
        split_v = v.split("_")
        result_len = result_len + int(split_v[3])
    return result_len
# cycle_gene_self  = open("cycle_gene_self.txt", "w")
#for line in blastin.readlines():
#    t = line.strip().split("\t")
#    if prev_seg != t[0] and prev_seg != "":
#        elen = prev_seg.split("_")[3]
#        if float(prev_len) / float(elen) > blast_ratio or prev_len > 2000:
#            blast_segs.add(prev_seg)
#        prev_seg = t[0]
#        prev_len = int(t[3])
#    else:
#        if float(t[2]) > 60:
#            prev_len = prev_len + int(t[3])
#        prev_seg = t[0]
#elen = int(prev_seg.split("_")[3])
#if float(prev_len) / float(elen) > blast_ratio or prev_len > 2000:
#    blast_segs.add(t[0])

for line in blastin.readlines():
    t = line.strip().split("\t")
#    if "EDGE_27711765_length_1685_cov_5" in prev_seg:
#        print(prev_len,"tdtd")
    if (prev_seg != t[0] and prev_seg != "") or (prev_ref != t[1] and prev_ref != ""):
#        if "EDGE_27711765_length_1685_cov_5" in prev_seg:
#            print(prev_len, elen,"xdxd")
        elen = prev_seg.split("_")[3]
        if float(prev_len) / float(elen) > blast_ratio or prev_len > 2000:
            blast_segs.add(prev_seg)
        prev_seg = t[0]
        prev_ref = t[1]
        prev_len = int(t[3])
    else:
        if float(t[2]) > blast_ratio*100:
#            if "EDGE_27711765_length_1685_cov_5" in prev_seg:
#                print(prev_len,float(t[2]),blast_ratio*100,"mmm")
            prev_len = prev_len + int(t[3])
        prev_seg = t[0]
        prev_ref = t[1]
elen = prev_seg.split("_")[3]
if float(prev_len) / float(elen) > blast_ratio or prev_len > 2000:
    blast_segs.add(t[0])
#print(blast_segs)
# gene_res = set()
# with open(gene_hit, 'r') as gene_lst:
#     for gene_r in gene_lst:
#         gene_res.add(gene_r.strip())


plasscore = {}
with open(plasscore_file, 'r') as ps:
    for s in ps:
        item = s.strip().split('\t')
        id_ = int(item[0].split('_')[1])
        plasscore[id_] = float(item[1])

genehit = []
# for id, i in enumerate(orderintmp):
    # genehit[id + 1] = False
with open(gene_hit, 'r') as gh:
    for s in gh:
        item = s.strip().split('\t')
        if len(s.strip()) == 0:
            continue
        genehit.append(s.strip())
        # id_ = int(item[0].split('_')[1])
        # genehit[id_] = True
#print(genehit)

record_dict = SeqIO.to_dict(SeqIO.parse(fain, "fasta"))
n_seq = Seq("N" * 40)
self_tag = False
cycle_tag = False
for idx, line in enumerate(orderin.readlines()):
    # print(idx, genehit[idx+1])
    if line.startswith("iter") or line.startswith("self"):
        if line.startswith("self"):
            self_tag = True
        if line.startswith("iter"):
            cycle_tag = True
        continue
    seq = ""
    tmp = line.strip().split("\t")
    if len(tmp) == 1 and self_tag:
        # if tmp[0][0:-1] not in blast_segs or plasscore[idx] < 0.7:
        #     continue
        if contains_gen(line, genehit) or plasscore[idx + 1] >= 0.7:
            for t in tmp:
                tmp_seq = record_dict[t[0:-1]].seq
                if t[-1] == '-':
                    tmp_seq = record_dict[t[0:-1]].seq.reverse_complement()
                seq = seq + tmp_seq
            #faout.write(">self-gene" + "".join(tmp) + "\n" + str(seq) + "\n")
            print("selfgene" + "".join(tmp))
            res_count.add('selfgene'+''.join(tmp))
            # cycle_gene_self.write(tmp+'\n')
        else:
            for t in tmp:
                tmp_seq = record_dict[t[0:-1]].seq
                if t[-1] == '-':
                    tmp_seq = record_dict[t[0:-1]].seq.reverse_complement()
                seq = seq + tmp_seq
            if "".join(tmp) not in in_faout:
                faout.write(">" + "".join(tmp) + "\n" + str(seq) + "\n")
                in_faout.append("".join(tmp))
            #res_count.add(''.join(tmp))
        continue

    if cycle_tag:
        # if tmp[0][0:-1] not in blast_segs or plasscore[idx] < 0.7:
        #     continue
        if contains_gen(line, genehit):
            print("cyclegene" + "".join(tmp))
            res_count.add('cyclegene'+''.join(tmp))
            # cycle_gene_self.write(tmp+'\n')
        if plasscore[idx + 1] >= 0.9:
            print("cyclescore" + "".join(tmp))
            res_count.add('cyclescore'+''.join(tmp))
            # cycle_gene_self.write(tmp+'\n')

    flags = False
    is_gene = False
    blast_len = 0
    all_len = 0
    if contains_gen(line, genehit):
        flags = True
        is_gene = True
        # if cycle_tag:
        #     print(">cycle-gene" + "".join(tmp))
    for t in tmp:
        elen = int(t.split("_")[3])
        all_len = all_len + elen
        if t[0:-1] in blast_segs:
            blast_len = blast_len + elen
    if (float(blast_len) / float(all_len)) > 0.2:
        flags = True
    if all_len < 1000:
        flags = False
    if not flags and plasscore[idx + 1] < 0.9 or all_len < 1000:
        continue
    for t in tmp:
        tmp_seq = record_dict[t[0:-1]].seq
        if t[-1] == '-':
            tmp_seq = record_dict[t[0:-1]].seq.reverse_complement()
        seq = seq + tmp_seq
    if contains_gen(line, genehit) and plasscore[idx + 1] >= 0.9:
        print("genescore" + "".join(tmp))
        if "".join(tmp) not in in_faout:
            faout.write(">" + "".join(tmp) + "\n" + str(seq) + "\n")
            in_faout.append("".join(tmp))
            #res_count.add('genescore'+''.join(tmp))
    else:
        if plasscore[idx + 1] >= 0.9:
            #print("score" + "".join(tmp))
            if "".join(tmp) not in in_faout:
                faout.write(">" + "".join(tmp) + "\n" + str(seq) + "\n")
                in_faout.append("".join(tmp))
        elif contains_gen(line, genehit):
            if "".join(tmp) not in in_faout:
                faout.write(">" + "".join(tmp) + "\n" + str(seq) + "\n")
                in_faout.append("".join(tmp))
            #print("gene" + "".join(tmp))
        if flags:
            #print(">ref" + "".join(tmp))
            if "".join(tmp) not in in_faout:
                faout.write(">" + "".join(tmp) + "\n" + str(seq) + "\n")
                in_faout.append("".join(tmp))

with open(sys.argv[8], 'w') as res:
    s_len=0
    for s in res_count:
        s_len = line_len(s)
        if s_len >= 10000:
            res.write(s+"\n")
