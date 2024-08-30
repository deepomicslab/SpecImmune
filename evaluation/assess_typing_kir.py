
import sys
import argparse
import pandas as pd

sys.path.insert(0, sys.path[0]+'/../scripts/')

from four_field_compare import assess_sim_module
from determine_gene import get_focus_gene



args = {}
args['i'] = 'KIR'
gene_list, interval_dict =  get_focus_gene(args)



outdir='/mnt/d/HLAPro_backup/Nanopore_optimize/data/sim_hap/reads3/'
resultdir='/mnt/d/HLAPro_backup/Nanopore_optimize/data/sim_hap/results3/'

dp=5
i=1

data = []
data2 = []
dp_list =  [5, 10, 15, 20, 25, 30, 35, 40]
# dp_list = [5, 10, 15, 20]
for dp in dp_list:
    data_dict = {}
    for i in range(1, 21):
        sample=f'KIR_dp{dp}_acc90_{i}'
        print (sample)
        true=f'{outdir}/{sample}/{sample}.KIR.hap.alleles.txt'
        infer= f'{resultdir}/{sample}/{sample}.KIR.final.type.result.txt' 
        spec_gene_accuracy_dict = assess_sim_module(true, infer, gene_list, args['i'])
        for gene in spec_gene_accuracy_dict:
            print (gene, spec_gene_accuracy_dict[gene])
            if gene not in data_dict:
                data_dict[gene] = [0, 0]
            data_dict[gene][0] += spec_gene_accuracy_dict[gene][1]
            data_dict[gene][1] += spec_gene_accuracy_dict[gene][2]
    
    accuracy_list = []
    for gene in gene_list:
        accuracy_list.append(data_dict[gene][0]/data_dict[gene][1])
        data2.append([dp*2, gene, data_dict[gene][0]/data_dict[gene][1]])
    print (accuracy_list)
    data.append(accuracy_list)

## transform the data to a pandas dataframe, and give row names
df = pd.DataFrame(data, columns=gene_list, index=dp_list)
print (df)
df.to_csv("kir_results/sim_results.csv", index=True)

df = pd.DataFrame(data2, columns=['Depth', 'Gene', 'Accuracy'])
df.to_csv("kir_results/sim_results2.csv", index=True)



