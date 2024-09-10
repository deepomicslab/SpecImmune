import os
import pandas as pd

def get_standard_diploid(pure_diplotype):
    diplotype_list = pure_diplotype.split("/")
    if len(diplotype_list) == 1:
        diplotype_list.append("NA")
    # print (pure_diplotype)
    return diplotype_list


def read_spec_result(spec_result):
    # check if the file exists
    spec_result_dict = {}
    spec_gene_depth = {}
    if not os.path.exists(spec_result):
        raise FileNotFoundError(f"{spec_result} does not exist")
    pure_diplotype = 'NA'
    with open(spec_result, 'r') as f:
        for line in f:
            field = line.strip().split("\t")
            # print (field)
            if field[1] == 'CYP2D6 Diplotype:':
                pure_diplotype = field[2]
            if field[1] == 'CYP2D6 Phenotype:':
                phenotype = field[2]
            if field[1] == 'CYP2D6 Activity_score:':
                activity = field[2]
            if field[0] != "#" and field[0] != "Locus":
                gene = field[0]
                allele = field[6]  #suballele
                read_num = int(field[4])
                spec_gene_depth[gene] = read_num
                if allele != "NA":
                    allele = allele.split(".")[0] # star allele
                    allele = "*" + allele.split('*')[1]

                if gene not in spec_result_dict:
                    spec_result_dict[gene] = []
                spec_result_dict[gene].append([allele])
    # print ("#", pure_diplotype, spec_result_dict)
    diplotype_list = get_standard_diploid(pure_diplotype)
    return diplotype_list, read_num, phenotype, activity
    

spec_dir = "/home/wangshuai/00.hla/long/experiments/cyp/cyp_results/spec_1k_all/"
result = "cyp_1k_all.csv"

data = []
for folder in os.listdir(spec_dir):
        ## check if the folder is a folder
        if not os.path.isdir(os.path.join(spec_dir, folder)):
            continue
        sample = folder
        spec_result = os.path.join(spec_dir, folder, f"{folder}.CYP.merge.type.result.txt")
        diplotype_list, read_num, phenotype, activity = read_spec_result(spec_result)
        data.append([sample,"CYP2D6", '1', diplotype_list[0], read_num, phenotype, activity])
        data.append([sample,"CYP2D6", '2', diplotype_list[1], read_num, phenotype, activity])

# transfer to dataframe
df = pd.DataFrame(data, columns=["Sample", "Locus", "Chromosome", "Genotype", "Reads_num", "Phenotype", "Activity_score"])
df.to_csv(result, index=False)