import sys
import re
inp = open(sys.argv[1])
inp2 = open(sys.argv[2])
outs = open(sys.argv[3],"w")
depth = sys.argv[4]
f_th = sys.argv[5]
blastin = open(sys.argv[6])
blast_ratio = float(sys.argv[7])
blast_segs = set()
relevate_blast_segs = set()
prev_seg = ""
prev_len = 0
relevate_edge_len = 2000
self_l_segs = set()
def get_len(edge):
    return int(edge.split("_")[3])
for line in blastin.readlines():
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
tmp = {}
for line in inp:
    vs = line.split("\t")
    a = re.split(":|,|;",vs[0])
    tmp[a[0]]=a[1:]
all_segs = {}
write_segs = set()
write_juncs = []
for line in inp2:
    vs = line.rstrip().split(" ")
    if vs[0] == "SEG":
        all_segs[vs[1]] = line
        # outs.write(line)
        continue
    if vs[1] == vs[3]:
        self_l_segs.add(vs[1])
        write_segs.add(all_segs[vs[1]])
        write_juncs.append(" ".join(vs)+"\n")
        continue
    if int(vs[-1]) < int(f_th) and (vs[1] not in blast_segs or vs[3] not in blast_segs):
        continue
    k = vs[1]
    kc = vs[1]+"'"
    if vs[2] == '-':
        k = vs[1]+"'"
        kc = vs[1]
    v = vs[3]
    vc = vs[3]+"'"
    if vs[4] == '-':
        v = vs[3]+"'"
        vc = vs[3]
    if k =="EDGE_5369_length_3828_cov_7.082378":
        print(tmp[k])
        print(v)
    if (k in tmp.keys() and v in tmp[k]):
        # print(v)
        if int(depth) > int(vs[-1]):
            vs[-1] = depth
        #write_juncs.append(" ".join(vs)+"\n")
        # outs.write(" ".join(vs)+"\n")
        tmp[k].remove(v)
    if (vc in tmp.keys() and kc in tmp[vc]):
        if int(depth) > int(vs[-1]):
            vs[-1] = depth
        #write_juncs.append(" ".join(vs)+"\n")
        # outs.write(" ".join(vs)+"\n")
        tmp[vc].remove(kc)
    if (vs[1] in blast_segs and vs[3] in blast_segs) or ((vs[1] in relevate_blast_segs or vs[3] in relevate_blast_segs) and float(vs[-1]) >= float(depth)/2 and get_len(vs[1]) <= relevate_edge_len and get_len(vs[3]) <= relevate_edge_len):
        if vs[1] not in blast_segs and (get_len(vs[1]) > relevate_edge_len or float(vs[-1]) <= float(depth)/2):
            continue
        if vs[3] not in blast_segs and (get_len(vs[3]) > relevate_edge_len or float(vs[-1]) <= float(depth)/2):
            continue
        if vs[1] in blast_segs and vs[2] in blast_segs:
            vs[-1] = str(int(vs[-1]) * 1.5)
        write_juncs.append(" ".join(vs)+"\n")
        write_segs.add(all_segs[vs[1]])
        write_segs.add(all_segs[vs[3]])
        # blast_segs.add(vs[1])
        # blast_segs.add(vs[3])
        relevate_blast_segs.add(vs[1])
        relevate_blast_segs.add(vs[2])
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
        if first == second:
            self_l_segs.add(first)
            write_segs.add(all_segs[first])
            write_juncs.append("JUNC {} {} {} {} {}\n".format(first,fdir,second,sdir,int(depth)/2))
            continue
        if (first in blast_segs and second in blast_segs) or (first in relevate_blast_segs or second in relevate_blast_segs and get_len(first) <= relevate_edge_len and get_len(second) <= relevate_edge_len):
            if first not in blast_segs and get_len(first) > relevate_edge_len:
                continue
            if second not in blast_segs and get_len(second) > relevate_edge_len:
                continue
            if first in blast_segs and second in blast_segs:
                write_juncs.append("JUNC {} {} {} {} {}\n".format(first,fdir,second,sdir,int(depth) * 2/3))
            else:
                write_juncs.append("JUNC {} {} {} {} {}\n".format(first,fdir,second,sdir,int(depth)/2))
            write_segs.add(all_segs[first])
            write_segs.add(all_segs[second])
            # blast_segs.add(first)
            # blast_segs.add(second)
            relevate_blast_segs.add(first)
            relevate_blast_segs.add(second)

writed_segs = []
for item in write_segs:
    writed_segs.append(item.split()[1])
    outs.write(item.strip()+" 0 1\n")
for item in blast_segs:
    if item not in writed_segs:
        outs.write(all_segs[item].strip()+" 0 1\n")
for item in write_juncs:
    t = item.split()
    if (t[1] in self_l_segs and t[3] not in self_l_segs) or (t[3] in self_l_segs and t[1] not in self_l_segs):
        continue
    outs.write(item)