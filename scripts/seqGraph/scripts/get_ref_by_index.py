import sys
import re
from Bio import SeqIO
fa = open(sys.argv[1])
fai = open(sys.argv[2])
ref_indx=open(sys.argv[3])
out_ref = open(sys.argv[4],"w")

fai_array = []
for line in fai.readlines():
    t = line.split("\t")
    fai_array.append(t[0])

record_dict = SeqIO.to_dict(SeqIO.parse(fa, "fasta"))

for line in ref_indx.readlines():
    t = re.split(r"\s+",line)
    print(t)
    ref_name = fai_array[int(t[1])-1]
    out_ref.write(">"+ref_name+"\n")
    out_ref.write(str(record_dict[ref_name].seq)+"\n")