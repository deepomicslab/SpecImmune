import os
import pandas as pd

def read_LD_values(indir,outfile):
    ## for each file in the indir
    ## read the LD values and store them in a df
    ## return the df
    data = []
    for file in os.listdir(indir):
        if file.endswith(".csv"):
            print(file)
            with open(indir + file, 'r') as f:
                i = 0
                for line in f:
                    if i == 1:
                        D = line.strip().replace('"', '')
                    if i == 2:
                        Wn = line.strip().replace('"', '')
                    i += 1
            field = file.split("_") 
            gene1 = field[0]
            gene2 = field[1]
            data.append([gene1, gene2, D, Wn])
            data.append([gene2, gene1, D, Wn])
            print (gene1, gene2, D, Wn)
            # break
    df = pd.DataFrame(data, columns = ['gene1', 'gene2', 'D', 'Wn'])
    df.to_csv(outfile, index = False)

# indir = "/mnt/d/HLAPro_backup/Nanopore_optimize/1kgp_analysis/kir_LD_result/"
# outfile = "kir_LD_values.csv"
# read_LD_values(indir,outfile)

# indir = "/mnt/d/HLAPro_backup/Nanopore_optimize/1kgp_analysis/hla_LD_result/"
# outfile = "hla_LD_values.csv"
# read_LD_values(indir,outfile)

indir = "/mnt/d/HLAPro_backup/Nanopore_optimize/1kgp_analysis/hla_kir_LD_result/"
outfile = "hla_kir_LD_values.csv"
read_LD_values(indir,outfile)
