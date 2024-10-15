import sys
from operator import itemgetter
import re

fai = open(sys.argv[1])  # fai
graph = open(sys.argv[2])  # graph
segf = open(sys.argv[4])  # 需要从新的segs
outs = sys.argv[3]  # 输出prefix
depth = float(sys.argv[5])  # samtools 计算而来
tmp = {}
segs = []
seg_g = {}
out_juncs = []
for line in fai:
    vs = line.split("\t")
    a = re.split(":|,|;", vs[0])
    tmp[a[0]] = a[1:]
ref_name_records = {}
for idx, line in enumerate(segf):
    l = line.strip().split('\t')
    t = re.split(r"[\+\-]", l[0])[:-1]
    print(idx, t)
    ref_name_records[idx] = l[1]
    line_segs = []
    out_juncs.append([])
    for v in t:
        line_segs.append(v)
        if v in seg_g.keys():
            if idx in seg_g[v]:
                continue
            seg_g[v].append(idx)
        else:
            seg_g[v] = [idx]
    segs.append(line_segs)
print(list(seg_g.values()))
out_segs = [[] for i in range(0, len(out_juncs))]


def in_segs(seg1, seg2=None):
    if seg2 is None:
        for s in segs:
            if seg1 in s:
                return 1
    else:
        res = []
        for i in range(len(segs)):
            s = segs[i]
            if seg1 in s and seg2 in s:
                res.append(i+1)
        return res
    return 0


for line in graph:
    vs = line.rstrip().split(" ")
    if vs[0] == "SEG":
        if in_segs(vs[1]):
            # print(vs, 'test', vs[1])
            for v in seg_g[vs[1]]:
                out_segs[v].append(line)
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
    # if vs[1] == "EDGE_72432_length_111_cov_3.857143" or vs[3] =="EDGE_72432_length_111_cov_3.857143":
        # print(in_segs(vs[1], vs[3]), "rrrrrrrrrr")
    # if vs ==[]:
        # print(in_segs(vs[1], vs[3]), 'dddddddd')

    if in_segs(vs[1], vs[3]):
        for it in in_segs(vs[1], vs[3]):
            out_juncs[it - 1].append(vs)
# for item in tmp.keys():
#     first = item
#     fdir = "+"
#     if item[-1] == "'":
#         first = item[0:-1]
#         fdir = "-"
#     for i in tmp[item]:
#         if i == "":
#             continue
#         second = i
#         sdir = "+"
#         if i[-1] == "'":
#             second = i[0:-1]
#             sdir = "-"
#         if in_segs(first, second):
#             vs = ["JUNC",first, fdir, second, sdir, str(int(depth) / 2)]
#             out_juncs[in_segs(vs[1], vs[3]) - 1].append(vs)
            # out_juncs.append(vs)
            # out_juncs("JUNC {} {} {} {} {}\n".format(first,fdir,second,sdir,int(depth)/2))
# s_juncs = sorted(out_juncs, key=itemgetter(2))
i = 0
# print(out_segs)
for js in out_segs:
    outs_f = open(outs + "_" + str(i) + "_" + ref_name_records[i] + ".second", "w")
    i = i + 1
    # print(i, js)
    for j in js:
        outs_f.write(j)
    outs_f.close()
i = 0
# print(out_juncs)
for js in out_juncs:
    print(i, js)
    outs_f = open(outs + "_" + str(i) + "_" + ref_name_records[i] + ".second", "a")
    i = i + 1
    for j in sorted(js):
        outs_f.write(" ".join(j) + "\n")
    outs_f.close()
