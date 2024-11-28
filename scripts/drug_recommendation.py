# write a script, use args to get allele type txt and drug db txt, and output the drug recommendation result
# Usage: python drug_recommendation.py -a allele_type.txt -d drug_db.txt
# Output: drug_recommendation.txt
# db txt header: Gene,Variant.Haplotypes,Genotype/Allele,Level.of.Evidence,Level.Modifiers（for metabolism）,Evidence_Score,Drug,Phenotype.Category,Phenotype.diseases,Phenotype.Genotype,Population,Allele Function,


import argparse
import pandas as pd
import numpy as np
import os

def get_args():
    parser = argparse.ArgumentParser(description='Drug recommendation based on allele type')
    parser.add_argument('-a', '--allele_type', help='allele type file', required=True)
    parser.add_argument('-d', '--drug_db', help='drug database file', required=True)
    parser.add_argument('-g', '--gene', help='gene_name', default='CYP', required=True)
    parser.add_argument('-o', '--output', help='output file', default='drug_recommendation.txt')
    return parser.parse_args()

def read_data(file):
    if not os.path.exists(file):
        print(f"Error: {file} not found")
        exit(1)
    # skip the first row
    data = pd.read_csv(file, sep='\t', skiprows=1)
    return data

def drug_recommendation(allele_type, drug_db):
    # Gene,Variant.Haplotypes,Genotype/Allele,Level.of.Evidence,Level.Modifiers（for metabolism）,Evidence_Score,Drug,Phenotype.Category,Phenotype.diseases,Phenotype.Genotype,Population,Allele Function,
    allele_type = read_data(allele_type)
    print(allele_type.head())
    drug_db = pd.read_csv(drug_db)
    # print(drug_db)
    # get all allele and their type
    allele_drug_dict = {}
    all_types=list(allele_type['Genotype'].unique())
    # remove nan
    all_types = [x for x in all_types if str(x) != 'nan']
    # split elelement in all_types with ";" in some cases to get all types
    print(len(all_types))
    print(all_types[:2])
    all_types = [x.split(';') for x in all_types]
    # merge all types in one list
    all_types = [item for sublist in all_types for item in sublist]
    # remove decimal part in the type
    all_types = [x.split('.')[0] for x in all_types]
    # unique all types
    all_types = set(all_types)
    print(list(all_types))
    # get drug for each type, if allele type is not in the drug db, return "No drug found"
    for allele in all_types:
        gene_name=allele.split('*')[0]
        # get rows with the same gene name
        gene_rows = drug_db[drug_db['Gene'] == gene_name]
        # search if haplotype contaning the Haplotypes
        for index, row in gene_rows.iterrows():
            allele_db=row['Genotype/Allele']
            star_allele="*"+allele.split('*')[1]

            if star_allele in allele_db.split('/'):
                if allele not in allele_drug_dict:
                    allele_drug_dict[allele] = []
                res=[row['Drug'], row['Level.of.Evidence'], row['Evidence_Score'], row['Level.Modifiers（for metabolism）']]
                if res not in allele_drug_dict[allele]:
                    allele_drug_dict[allele].append(res)
        
    print(allele_drug_dict)
    res_f=open(f'{args.output}', 'w')
    for allele in allele_drug_dict:
        print(f"Allele: {allele}")
        res_f.write(f"Allele: {allele}\n")
        for drug in allele_drug_dict[allele]:
            print(f"\tDrug: {drug[0]}, Level of Evidence: {drug[1]}, Evidence Score: {drug[2]}, Level Modifiers: {drug[3]}")
            res_f.write(f"\tDrug: {drug[0]}, Level of Evidence: {drug[1]}, Evidence Score: {drug[2]}, Level Modifiers: {drug[3]}\n")
    res_f.close()
        


if __name__ == '__main__':
    args = get_args()
    drug_recommendation(args.allele_type, args.drug_db)

