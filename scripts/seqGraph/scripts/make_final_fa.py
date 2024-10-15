from Bio import SeqIO
import re
import sys
from Bio.Seq import Seq

fain = sys.argv[1]
orderin = open(sys.argv[2])
faout= open(sys.argv[3], "w")
prefix = sys.argv[4]

record_dict = SeqIO.to_dict(SeqIO.parse(fain, "fasta"))
n_seq = Seq("N"*10)

count = 0
for line in orderin.readlines():
    if "all" in line:
        continue
    seq = ""
    tmp = re.split("\t+", line.strip("\n"))
    for t in tmp:
        if t == '':
            continue
        t = t.replace("ref","")
        tmp_seq = record_dict[t[0:-1]].seq
        if t[-1] == '-':
            tmp_seq = record_dict[t[0:-1]].seq.reverse_complement()
        seq = seq + tmp_seq
    count += 1
    faout.write(">" + prefix + "_phage_" + str(count) + "\n")
    faout.write(str(seq) + "\n")

orderin.close()
faout.close()
