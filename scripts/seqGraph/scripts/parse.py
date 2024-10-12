import sys
import re

inp = open(sys.argv[1])
inp2 = open(sys.argv[2])
outs = open(sys.argv[3], "w")
depth = sys.argv[4]
f_th = sys.argv[5]
gene_file = sys.argv[6]
score_file = sys.argv[7]
blast_in = open(sys.argv[8])
blast_ratio = float(sys.argv[9])
to_remove_score_threhold = 0.2
tmp = {}
gene_res = {}
scores = {}

blast_segs = set()
relevate_blast_segs = set()
prev_seg = ""
prev_len = 0
for line in blast_in.readlines():
    t = line.strip().split("\t")
    if prev_seg != t[0] and prev_seg != "":
        elen = prev_seg.split("_")[3]
        if float(prev_len)/float(elen) > blast_ratio or prev_len > 2000:
            blast_segs.add(prev_seg)
        prev_seg = t[0]
        prev_len = int(t[3])
    else:
        if float(t[2]) > 60:
            prev_len = prev_len + int(t[3])
        prev_seg = t[0]
elen = prev_seg.split("_")[3]
if float(prev_len)/float(elen) > blast_ratio or prev_len > 2000:
    blast_segs.add(t[0])
# print(blast_segs)
if len(sys.argv) > 5:
    with open(gene_file, 'r') as gene_lst:
        for gene_r in gene_lst:
            gene_res[gene_r.strip()] = '1'

    with open(score_file, 'r') as score_lst:
        for score_r in score_lst:
            line_lst = score_r.strip().split("\t")
            scores[line_lst[0]] = str(line_lst[1])

for line in inp:
    vs = line.split("\t")
    a = re.split(":|,|;", vs[0])
    tmp[a[0]] = a[1:]
# print(tmp)
all_segs = {}
write_segs = set()
write_juncs = []
for line in inp2:
    vs = line.rstrip().split(" ")

    if vs[0] == "SEG":
        all_segs[vs[1]] = line
        # outs.write(line)
        continue
    if int(vs[-1]) < int(f_th):
        continue
    # print(vs)
    left_score = float(scores[vs[1]])
    right_score = float(scores[vs[3]])
    left_node_len = int(vs[1].split("_")[3])
    right_node_len = int(vs[3].split("_")[3])
    if (left_score < 0.2 and left_node_len > 10000) or (right_score < 0.2 and right_node_len > 10000):
        # print(vs, left_score, right_score)
        continue
    k = vs[1]
    kc = vs[1] + "'"
    if vs[2] == '-':
        k = vs[1] + "'"
        kc = vs[1]
    v = vs[3]
    vc = vs[3] + "'"
    if vs[4] == '-':
        v = vs[3] + "'"
        vc = vs[3]
    # if k =="EDGE_5369_length_3828_cov_7.082378":
    #     print(tmp[k])
    #     print(v)
    # print("k:", k)
    # print("kc", kc)
    # print("v", v)
    # print("vc", vc)
    is_add_junc = False
    if (k in tmp.keys() and v in tmp[k]):
        # print(v)
        if int(depth) > int(vs[-1]):
            vs[-1] = depth
        # write_juncs.append(" ".join(vs)+"\n")
        # outs.write(" ".join(vs)+"\n")
        # write_juncs.append(" ".join(vs)+"\n")
        tmp[k].remove(v)
    if (vc in tmp.keys() and kc in tmp[vc]):
        if int(depth) > int(vs[-1]):
            vs[-1] = depth
        # write_juncs.append(" ".join(vs)+"\n")
        # outs.write(" ".join(vs)+"\n")
        # write_juncs.append(" ".join(vs)+"\n")
        tmp[vc].remove(kc)
    if vs[1] in blast_segs or vs[3] in blast_segs or vs[1] in gene_res or vs[3] in gene_res or vs[1] == vs[3] or (float(vs[-1]) >= float(depth)/4 and (vs[1] in relevate_blast_segs or vs[3] in relevate_blast_segs)):
        write_juncs.append(" ".join(vs) + "\n")
    # print(all_segs[vs[1]].strip()+" "+(gene_res[vs[1]] if vs[1] in gene_res else '0')+" "+(scores[vs[1]] if vs[1] in scores else '0')+"\n")
        write_segs.add(all_segs[vs[1]].strip() + " " + (gene_res[vs[1]] if vs[1] in gene_res else '0') + (
            " " + scores[vs[1]] if vs[1] in scores else '0') + "\n")
        write_segs.add(all_segs[vs[3]].strip() + " " + (gene_res[vs[3]] if vs[3] in gene_res else '0') + " " + (
            scores[vs[3]] if vs[3] in scores else '0') + "\n")
        relevate_blast_segs.add(vs[1])
        relevate_blast_segs.add(vs[3])
# print(tmp['EDGE_5369_length_3828_cov_7.082378'])

for item in tmp.keys():
    first = item
    fdir = "+"
    if item[-1] == "'":
        first = item[0:-1]
        fdir = "-"
    for i in tmp[item]:
        if i == "":
            continue
        second = i
        sdir = "+"
        if i[-1] == "'":
            second = i[0:-1]
            sdir = "-"

        left_score = float(scores[first])
        right_score = float(scores[second])
        left_node_len = int(first.split("_")[3])
        right_node_len = int(second.split("_")[3])
        if ((left_score < 0.2 and left_node_len > 10000) or (right_score < 0.2 and right_node_len > 10000)) and \
            (first not in blast_segs and second not in blast_segs and first not in gene_res and second not in gene_res and first != second) and (first not in relevate_blast_segs or second not in relevate_blast_segs):
            # print(first, left_score, second, right_score)
            continue
        # if first in blast_segs or second in blast_segs or first in gene_res or second in gene_res:
        write_juncs.append("JUNC {} {} {} {} {}\n".format(first, fdir, second, sdir, int(depth) / 2))
        write_segs.add(all_segs[first].strip() + " " + (gene_res[first] if first in gene_res else '0') + " " + (
            scores[first] if first in scores else '0') + "\n")
        write_segs.add(all_segs[second].strip() + " " + (gene_res[second] if second in gene_res else '0') + " " + (
            scores[second] if second in scores else '0') + "\n")
        relevate_blast_segs.add(first)
        relevate_blast_segs.add(second)
for item in write_segs:
    # print(item)
    outs.write(item)
for item in write_juncs:
    outs.write(item)
