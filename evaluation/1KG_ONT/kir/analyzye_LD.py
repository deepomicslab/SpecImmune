import os
import pandas as pd

def count_line_number_of_file(file):
    with open(file, 'r') as f:
        i = 0
        for line in f:
            i += 1
    return i 

def read_LD_values(raw_dir, indir,outfile):
    ## for each file in the indir
    ## read the LD values and store them in a df
    ## return the df
    data = []
    gene_set = set()
    for file in os.listdir(indir):
        if file.endswith(".csv"):
            # print(file)
            with open(indir + file, 'r') as f:
                i = 0
                for line in f:
                    if i == 1:
                        D = line.strip().replace('"', '')
                    if i == 2:
                        Wn = line.strip().replace('"', '')
                    if i == 3:
                        Wab = line.strip().replace('"', '')
                    if i == 4:
                        Wba = line.strip().replace('"', '')
                    if i == 5:
                        hap_num = line.strip().replace('"', '')
                    i += 1
            sample_num = count_line_number_of_file(raw_dir + file) -1
            if Wn == "Inf" or Wn == "NaN":
                continue
            field = file.split("_") 
            gene1 = field[0]
            gene2 = field[1]
            min_w = min(float(Wab), float(Wba))
            data.append([gene1, gene2, D, Wn, hap_num, sample_num, Wab, min_w])
            data.append([gene2, gene1, D, Wn, hap_num, sample_num, Wba, min_w])
            gene_set.add(gene1)
            gene_set.add(gene2)

            # print (gene1, gene2, D, Wn)
            # break
    for gene in gene_set:
        data.append([gene, gene, 1, 1, 1000, 1000, 1, 1])
    df = pd.DataFrame(data, columns = ['gene1', 'gene2', 'D', 'Wn', 'hap_num','sample_num','Wab', 'min_w'])
    df.to_csv(outfile, index = False)

# indir = "/mnt/d/HLAPro_backup/Nanopore_optimize/1kgp_analysis/kir_LD_result/"
# outfile = "kir_LD_values.csv"
# read_LD_values(indir,outfile)

# indir = "/mnt/d/HLAPro_backup/Nanopore_optimize/1kgp_analysis/hla_LD_result/"
# outfile = "hla_LD_values.csv"
# read_LD_values(indir,outfile)

# raw_dir = "/mnt/d/HLAPro_backup/Nanopore_optimize/1kgp_analysis/hla_kir_LD/"
# indir = "/mnt/d/HLAPro_backup/Nanopore_optimize/1kgp_analysis/hla_kir_LD_result/"
# outfile = "hla_kir_LD_values.csv"
# read_LD_values(raw_dir, indir,outfile)

# raw_dir = "/mnt/d/HLAPro_backup/Nanopore_optimize/1kgp_analysis/hla_kir_cyp_LD/"
# indir = "/mnt/d/HLAPro_backup/Nanopore_optimize/1kgp_analysis/hla_kir_cyp_LD_result/"
# outfile = "hla_kir_cyp_LD_values.csv"
# read_LD_values(raw_dir, indir,outfile)


raw_dir = "/mnt/d/HLAPro_backup/Nanopore_optimize/1kgp_analysis/hla_kir_cyp_vdj_LD/"
indir = "/mnt/d/HLAPro_backup/Nanopore_optimize/1kgp_analysis/hla_kir_cyp_vdj_LD_result/"
outfile = "hla_kir_cyp_vdj_LD_values.csv"
read_LD_values(raw_dir, indir,outfile)
