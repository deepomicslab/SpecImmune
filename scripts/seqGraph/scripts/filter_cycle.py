import re
import os
import sys
res = set()
with open(sys.argv[1]) as r:
    for line in r.readlines():
        line_len = 0
        splited=re.split(r'[+-]',line.strip())
        for v in splited:
            if v == "" or v ==" ":
                continue
            line_len = line_len + int(v.split('_')[3])
        if line_len >= 10000:
        #print(line,"xxx\n")
            liner = line.replace("cycle","").replace("score","").replace("self","").replace("gene","").replace("ref","")
            res.add(liner.strip("\n"))
for item in res:
    print(item.replace("+","+\t").replace("-","-\t"))
