from Bio import SeqIO
import sys
from Bio.Seq import Seq
fain = sys.argv[1]
orderin = open(sys.argv[2])
faout= open(sys.argv[3],"w")
blastin = open(sys.argv[4])
blast_ratio = float(sys.argv[5])
blast_segs = set()
prev_seg = ""
prev_len = 0
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

record_dict = SeqIO.to_dict(SeqIO.parse(fain, "fasta"))
n_seq = Seq("N"*40)
print(blast_segs)
for line in orderin.readlines():
    if line.startswith("iter") or line.startswith("self"):
        continue
    seq = ""
    tmp = line.strip().split("\t")
    if len(tmp) == 1:
        if tmp[0][0:-1] not in blast_segs:
            continue
    flags = False
    blast_len = 0
    all_len = 0
    for t in tmp:
        elen = int(t[0:-1].split("_")[3])
        all_len = all_len + elen
        if t[0:-1] in blast_segs:
            blast_len = blast_len + elen
    if (float(blast_len) / float(all_len)) > 0.1:
        flags = True
    if flags and all_len < 1000:
        flags = False
    if not flags:
        continue
    for t in tmp:
        tmp_seq = record_dict[t[0:-1]].seq
        if t[-1] == '-':
            tmp_seq = record_dict[t[0:-1]].seq.reverse_complement()
        seq = seq + tmp_seq
    faout.write(">"+"".join(tmp)+"\n"+str(seq)+"\n")
