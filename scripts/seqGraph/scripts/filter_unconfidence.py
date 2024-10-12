import sys
import re
unconfidence_file = open(sys.argv[1])
confidence_file = open(sys.argv[2])
gene_file = open(sys.argv[3])

def get_len_arr(array):
    result = []
    for item in array:
        result.append(int(item.split("_")[3]))
    result.sort(reverse = True)
    return result

def line_gene_segs_num(line_arr, genes):
    counts = 0
    for item in line_arr:
        if item in genes:
            counts = counts + 1
    return counts

# def compare_two(line1,line2):
#     line1_arr = line1.strip().replace("+","").replace("-","").split("\t")
#     line1_lens = get_len_arr(line1_arr)
#     line2_arr = line2.strip().replace("+","").replace("-","").split("\t")
#     line2_lens = get_len_arr(line2_arr)
#     same_items = []
#     for item1 in line1_lens:
#         if item1 in line2_lens:
#             same_items.append(item1)
#     if sum(same_items)/sum(array1) >= 0.75 and sum(same_items)/sum(array2) >= 0.75 :
#         return 1
#     else:
#         return 0

def compare_two(line1,line2,genes):
    line1_arr = line1.strip().replace("+","").replace("-","").split("\t")
    line1_lens = get_len_arr(line1_arr)
    line2_arr = line2.strip().replace("+","").replace("-","").split("\t")
    line2_lens = get_len_arr(line2_arr)
    same_items = []
    for item1 in line1_lens:
        if item1 in line2_lens:
            same_items.append(item1)
    if (sum(same_items)/sum(line1_lens) >= 0.75 and sum(same_items)/sum(line2_lens) >= 0.75) or sum(same_items)/sum(line1_lens) >= 0.85 or sum(same_items)/sum(line2_lens) >= 0.85 :
        if line_gene_segs_num(line1_arr, genes) > line_gene_segs_num(line2_arr, genes):
            return 0
        elif line_gene_segs_num(line1_arr, genes) < line_gene_segs_num(line2_arr, genes):
            return 1
        else:
            if sum(line1_lens) >= sum(line2_lens):
                return 0
            else:
                return 1
    else:
        return 2
confidence_segs = []
# confidence
for line in confidence_file.readlines():
    line_arr = line.strip().replace("+","").replace("-","").split("\t")
    for item in line_arr:
        confidence_segs.append(item)
filtered_lines = []
# remove lines that all segs in confidence
for line in unconfidence_file.readlines():
    flag = False
    line_arr = line.strip().replace("+","\t").replace("-","\t").split("\t")
    for item in line_arr:
        itemr = item.replace("gene","").replace("score","").replace("cycle","")
        if itemr in confidence_segs:
            flag = True
            break
    if not flag:
        filtered_lines.append(line.strip())

# remove duplicate
# dsu = DSU(len(filtered_lines))
genes = []
for line in gene_file.readlines():
    genes.append(line.strip())
def remove_dup(lines):
    skip_idx = []
    final_lines = set()
    for i in range(0,len(lines)):
        need_add = True
        if i in skip_idx:
            continue
        line = lines[i]
        for j in range(i+1,len(lines)):
            if j in skip_idx:
                continue
            line2 = lines[j]
            tag = compare_two(line,line2,genes)
            if tag == 0:
                #final_lines.add(line)
                skip_idx.append(j)
                need_add = True
            elif tag == 1:
                #final_lines.add(line2)
                need_add = False
                break
        if need_add:
            final_lines.add(line)
    return list(final_lines)
prev_len = 0
while True:
    if prev_len == len(filtered_lines):
        break
    else:
        prev_len = len(filtered_lines)
        filtered_lines = remove_dup(filtered_lines)
        #print(filtered_lines)

for line in filtered_lines:
    line_len = 0
    line_re = re.split(r"[+-]",line.strip())
    for item in line_re:
        if item == "":
            continue
        line_len = line_len + int(item.split("_")[3])
    if line_len >= 10000:
        print(line)

# for line in filtered_lines:
#     line_arr = line.strip().replace("+","").replace("-","").split("\t")
#     line_lens = get_len_arr(line_arr)
#     for line2 in filtered_lines:
#         line2_arr = line2.strip().replace("+","").replace("-","").split("\t")
#         line2_lens = get_len_arr(line2_arr)
