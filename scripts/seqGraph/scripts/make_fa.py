from Bio import SeqIO
import sys
import re
def print_help():
    print("Make fasta from matching result.\nusage:\n\tpython make_fa.py row.fa result.txt out.fa\n"
          "\trow.fa: row assembly fasta\n\tresult.txt: matching result from program matching\n\tout.fa: output fasta")
fain = sys.argv[1]
if fain == "help":
    print_help()
    exit(0)
oin = sys.argv[2]
orderin = open(oin)
faout_f= sys.argv[3]
faout = open(faout_f,"w")
record_dict = SeqIO.to_dict(SeqIO.parse(fain, "fasta"))
h_appends = []

for line in orderin.readlines():
    if line.startswith("iter") or line.startswith("self"):
        continue
    seq = ""
    tmp = re.split("\t",line.strip())
    for t in tmp:
        tmp_seq = str(record_dict[t[0:-1]].seq)
        if t[-1] == '-':
            tmp_seq = str(record_dict[t[0:-1]].seq)[::-1]
        seq = seq + tmp_seq
        h_appends.append(t[0:-1])
        del record_dict[t[0:-1]]
    faout.write(">"+"_".join(tmp)+"\n")
    faout.write(seq+"\n")
for l in record_dict.keys():
    faout.write(">"+l+"\n")
    faout.write(str(record_dict[l].seq)+"\n")