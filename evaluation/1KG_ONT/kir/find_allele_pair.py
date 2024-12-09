
def get_gene_list(file):
    with open(file, 'r') as f:
        for line in f:
            field = line.strip().split(",")
            sample = field[0]
            sample_type = {}
            for allele in  field[1:]:
                gene = allele.split("*")[0]
                if gene not in sample_type:
                    sample_type[gene] = []
            break
    return list(sample_type.keys())

def read_csv(file, focus_gene1, focus_gene2):
    allele_pair_count = {}
    with open(file, 'r') as f:
        for line in f:
            field = line.strip().split(",")
            sample = field[0]
            sample_type = {}
            for allele in  field[1:]:
                gene = allele.split("*")[0]
                if gene not in sample_type:
                    sample_type[gene] = []
                sample_type[gene].append(allele)
            if focus_gene1 in sample_type and focus_gene2 in sample_type:
                for allele1 in sample_type[focus_gene1]:
                    for allele2 in sample_type[focus_gene2]:
                        pair = '&'.join(sorted([allele1, allele2]))
                        if pair not in allele_pair_count:
                            allele_pair_count[pair] = 0
                        allele_pair_count[pair] += 1
    ### sort the allele pair by the count in reverse order
    allele_pair_count = dict(sorted(allele_pair_count.items(), key=lambda item: item[1], reverse=True))
    ## print top 10 allele pair
    i = 0
    for pair, count in allele_pair_count.items():
        if count > 20:
            print (pair, count)
        i += 1
        if i == 10:
            break

file = "../chromoMap/sample_all_loci_type2.csv"
gene_list = get_gene_list(file)

# for i in range(len(gene_list)):
#     for j in range(i+1, len(gene_list)):
#         if gene_list[i][:2] == 'IG' or gene_list[j][:2] == 'IG':
#             continue
#         if gene_list[i][:2] == 'TR' or gene_list[j][:2] == 'TR':
#             continue
#         if  (gene_list[i][:3] == 'HLA' and gene_list[j][:3] == 'KIR'):
#             read_csv(file, gene_list[i], gene_list[j])
#         elif (gene_list[i][:3] == 'KIR' and gene_list[j][:3] == 'HLA'):
#             read_csv(file, gene_list[j], gene_list[i])

read_csv(file, 'HLA-A', 'KIR2DL1')
read_csv(file, 'HLA-DPB1', 'KIR2DL1')
