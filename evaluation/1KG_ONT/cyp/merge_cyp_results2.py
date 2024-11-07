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
                spec_result_dict[gene].append([field])
    # print ("#", pure_diplotype, spec_result_dict)
    diplotype_list = get_standard_diploid(pure_diplotype)
    return diplotype_list, read_num, phenotype, activity, spec_result_dict, spec_gene_depth
    

spec_dir = "/home/wangshuai/00.hla/long/experiments/cyp/cyp_results/spec_1k_all2/"
result = "cyp_1k_all2.csv"

data = []
for folder in os.listdir(spec_dir):
    ## check if the folder is a folder
    if not os.path.isdir(os.path.join(spec_dir, folder)):
        continue
    sample = folder
    spec_result = os.path.join(spec_dir, folder, f"{folder}.CYP.merge.type.result.txt")
    diplotype_list, read_num, phenotype, activity, spec_result_dict, spec_gene_depth = read_spec_result(spec_result)
    data.append(["CYP2D6", '1', diplotype_list[0], '-', read_num, '-', '-', sample])
    data.append(["CYP2D6", '2', diplotype_list[1], '-', read_num, '-', '-', sample])
    for gene in spec_result_dict:
        if gene == 'CYP2D6':
            continue
        for i in range(2):
            field = spec_result_dict[gene][i]
            array = field[2].split(";")
            for j in range(len(array)):
                if array[j] == "NA":
                    continue
                array[j] = array[j].split(".")[0]  ## only get star allele
            field[2] = ";".join(array)
            field[6] = field[6].split(".")[0]  ## only get star allele
            data.append([gene, i+1, field[2], field[3], field[4], field[5], field[6], sample])



# transfer to dataframe
df = pd.DataFrame(data, columns=["Locus","Chromosome","Genotype","Match_info","Reads_num","Step1_type","One_guess","Sample"])
df.to_csv(result, index=False)