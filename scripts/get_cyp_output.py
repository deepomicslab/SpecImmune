import json
import os
import sys

from cyp_phenotyper import phenotyper

cnv_alleles = ['*68', '*61', '*63', '*4.013', '*36', '*83', '*10', '*17', '*13', '*79','*80', '*78', '*67', '*66', '*76','*5']

## read a jason file into a dictionary
def read_pangu_result(pangu_result):
    # check if the file exists
    if not os.path.exists(pangu_result):
        raise FileNotFoundError(f"{pangu_result} does not exist")
    with open(pangu_result, 'r') as f:
        pangu_result = json.load(f)

    diplotype = pangu_result[0]['diplotype']
    # print (diplotype)
    # print (pangu_result[0]['haplotypes'])
    pangu_alleles = set()
    for hap_dict in pangu_result[0]['haplotypes']:
        for allele_dict in hap_dict['alleles']:
            pangu_alleles.add(allele_dict['call'])
    # print (alleles)
        
    return pangu_result, pangu_alleles


## read speclong result
def read_spec_result(spec_result_file):
    # check if the file exists
    if not os.path.exists(spec_result_file):
        raise FileNotFoundError(f"{spec_result_file} does not exist")
    with open(spec_result_file, 'r') as f:
        spec_result = []
        for line in f:
            field = line.strip().split()
            spec_result.append(field)
    return spec_result

def check_cnv(pangu_alleles, cnv_alleles):
    for allele in pangu_alleles:
        if allele in cnv_alleles:
            return True
    return False

def merge_result(pangu_result, pangu_alleles, spec_result, read_cutoff=10):
    output_result = spec_result.copy()
    detailed_diplotype = 'NA'
    diplotype = pangu_result[0]['diplotype'].split()[1]
    # print (diplotype)
    if len(pangu_alleles) > 2:
        print ("cnv detected, use pangu result")
    elif check_cnv(pangu_alleles, cnv_alleles):
        print ("cnv detected, use pangu result")
        
    else:  # use spec result
        spec_alleles = []
        copynumber = pangu_result[0]['copynumber']
        for i in [2, 3]:
            # print (spec_result[i][4])
            if int(spec_result[i][4]) < read_cutoff and copynumber == 2: # use step 1 result
                spec_alleles.append(spec_result[i][5])
            else:
                if spec_result[i] != 'NA':
                    spec_alleles.append(spec_result[i][2].split(';')[0])
        spec_alleles, star_alleles = refine_spec_alleles(spec_alleles)
            
        if copynumber <= 2:
            diplotype = "/".join(star_alleles)
            detailed_diplotype = "/".join(spec_alleles)
        else:
            diplotype, detailed_diplotype = update_diploid(star_alleles, spec_alleles, diplotype)
    
    diplotype_array = diplotype.split("/") 
    if len(diplotype_array) == 2:
        phenotype = phenotyper("cyp2d6", diplotype_array[0], diplotype_array[1])
    elif len(diplotype_array) == 1:
        phenotype = phenotyper("cyp2d6", diplotype_array[0], 'NA')
    else:
        print ("WARNING: diplotype is not correct")
        phenotype = 'NA'
    # print(phenotype)
    print (diplotype, detailed_diplotype, phenotype)
    return diplotype, detailed_diplotype, phenotype

def refine_spec_alleles(spec_alleles):
    if len(spec_alleles) == 1:
        spec_alleles = spec_alleles + spec_alleles
    star_alleles  =  []
    for i in range(len(spec_alleles)):
        spec_alleles[i] = "*" + spec_alleles[i].split("*")[1]
        star_alleles.append("*" + spec_alleles[i].split("*")[1].split('.')[0])
    return spec_alleles, star_alleles

def extract_pangu_diploid(diplotype):
    field = diplotype.split("/")
    suffix = []
    star = []
    for i in range(len(field)):
        fie = field[i].split("x")
        if len(fie) == 2:
            suffix.append("x"+fie[1])
        else:
            suffix.append('')
        star.append(fie[0])
    return suffix, star

def update_diploid(star_alleles, spec_alleles, diplotype):
    detailed_diplotype = 'NA'
    suffix, star = extract_pangu_diploid(diplotype)
    if star_alleles[0] == star[0] and star_alleles[1] == star[1]:
        pass
    elif star_alleles[0] == star[1] and star_alleles[1] == star[0]:
        ## reverse a list
        suffix = suffix[::-1]
    elif star_alleles[0] == star[0] and star_alleles[1] != star[1]:
        pass
    elif star_alleles[0] != star[0] and star_alleles[1] == star[1]:
        pass
    elif star_alleles[0] == star[1] and star_alleles[1] != star[0]:
        suffix = suffix[::-1]
    elif star_alleles[0] != star[1] and star_alleles[1] == star[0]:
        suffix = suffix[::-1]
    else:
        return diplotype, detailed_diplotype
    diplotype = star_alleles[0] + "" + suffix[0] + "/" + star_alleles[1] + "" + suffix[1]
    detailed_diplotype = spec_alleles[0]+ "" + suffix[0] + "/" + spec_alleles[1]+ "" + suffix[1]
    return diplotype, detailed_diplotype


def output(spec_result, diplotype, detailed_diplotype, phenotype, merge_result_file):
    f = open(merge_result_file, 'w')
    for i in range(1):
        print ('\t'.join(spec_result[i]), file = f)
    print ("#\tDiplotype\t"+diplotype, file = f)
    print (f"#\tPhenotype:\t{phenotype[1]}", file = f)
    print (f"#\tActivity_score:\t{phenotype[0]}", file = f)
    print ("#\tDetailed_diplotype\t"+detailed_diplotype, file = f)
    for i in range(1, len(spec_result)):
        print ('\t'.join(spec_result[i]), file = f)
    f.close()


### to do
# refine amplicon output







if __name__ == "__main__":  

    # prefix = "/mnt/d/HLAPro_backup/Nanopore_optimize/cyp_results/HG00436_1/HG00436_1"
    prefix = sys.argv[1]
    pangu_result = f"{prefix}_report.json"
    spec_result_file = f"{prefix}.CYP.final.type.result.txt"
    merge_result_file = f"{prefix}.CYP.merge.type.result.txt"
    pangu_result, pangu_alleles = read_pangu_result(pangu_result)
    spec_result = read_spec_result(spec_result_file)
    diplotype, detailed_diplotype, phenotype = merge_result(pangu_result, pangu_alleles, spec_result)
    output(spec_result, diplotype, detailed_diplotype, phenotype, merge_result_file)