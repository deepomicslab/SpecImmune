import json
import os
import sys
from collections import defaultdict

from cyp_phenotyper import phenotyper

cnv_alleles = ['*68', '*61', '*63', '*4.013', '*36', '*83', '*10', '*17', '*13', '*79','*80', '*78', '*67', '*66', '*76','*5']

## read a jason file into a dictionary
def read_pangu_result(pangu_result):
    # check if the file exists
    if not os.path.exists(pangu_result):
        print (f"WARNING: {pangu_result} does not exist")
        return {}, set()
    with open(pangu_result, 'r') as f:
        pangu_result = json.load(f)

    # diplotype = pangu_result[0]['diplotype']
    # print (diplotype)
    # print (pangu_result[0]['haplotypes'])
    pangu_alleles = set()
    if len(pangu_result) >= 1:
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
    spec_result = []
    all_spec_result = []
    with open(spec_result_file, 'r') as f:
        i = 0
        for line in f:
            field = line.strip().split()
            if i <= 1:
                spec_result.append(field)
            if field[0] == 'CYP2D6':
                spec_result.append(field)
            all_spec_result.append(field)
            i += 1
    return spec_result, all_spec_result

def check_cnv(pangu_alleles, cnv_alleles):
    for allele in pangu_alleles:
        if allele in cnv_alleles:
            return True
    return False

def check_36_10(diplotype):
    ## if it only has one haplotype, and it contains 36+10, then use spec result
    if diplotype == '*36+*10' or diplotype == '*36x2+*10' or diplotype == '*36x3+*10' or diplotype == '*36x4+*10':
        return True
    return False

def merge_result(pangu_result, pangu_alleles, spec_result, over2hap, tool='spec', read_cutoff=10):
    output_result = spec_result.copy()
    detailed_diplotype = 'NA'
    if len(pangu_result) >= 1:
        diplotype = pangu_result[0]['diplotype'].split()[1]
    # print (diplotype, pangu_alleles)
    if len(pangu_result) >= 1 and len(pangu_alleles) > 2:
        print ("cnv detected, use pangu result")
    elif len(pangu_result) >= 1 and  check_cnv(pangu_alleles, cnv_alleles) and not check_36_10(diplotype):
        print ("cnv detected, use pangu result")
    elif len(pangu_result) >= 1 and over2hap:
        print ("over 2 hap detected, use pangu result")
    elif len(pangu_result) >= 1 and seq_data == 'pacbio-hifi' and seq_tech != 'amplicon':
        print ("wgs hifi data, use pangu result")
    else:  # use spec result
        print ("no cnv detected, use spec result")
        
        if len(pangu_result) >= 1:
            copynumber = pangu_result[0]['copynumber']
        else:
            copynumber = 2
        if tool == 'spec':
            
            spec_alleles = []
            for i in [2, 3]:
                # print (spec_result[i][4])
                if int(spec_result[i][4]) < read_cutoff and copynumber == 2: # use step 1 result
                    spec_alleles.append(spec_result[i][5].split(';')[0])
                elif spec_result[i][2] == 'NA':
                    spec_alleles.append(spec_result[i][5].split(';')[0])
                else:
                    # if spec_result[i][2] != 'NA':
                    spec_alleles.append(spec_result[i][2].split(';')[0])
            # print (spec_alleles)
            spec_alleles, star_alleles = refine_spec_alleles(spec_alleles)
        else:  ## for stargazer
            spec_alleles, star_alleles = spec_result, spec_result
        print (spec_alleles, star_alleles)    
        if copynumber <= 2:
            diplotype = "/".join(star_alleles)
            detailed_diplotype = "/".join(spec_alleles)
        else:
            diplotype, detailed_diplotype = update_diploid(star_alleles, spec_alleles, diplotype)
    
    # print(phenotype)
    print (diplotype, detailed_diplotype)
    return diplotype, detailed_diplotype

def get_phenotype(diplotype):
    diplotype_array = diplotype.split("/") 
    if len(diplotype_array) == 2:
        phenotype = phenotyper("cyp2d6", diplotype_array[0], diplotype_array[1])
    elif len(diplotype_array) == 1:
        phenotype = phenotyper("cyp2d6", diplotype_array[0], 'NA')
    else:
        print ("WARNING: diplotype is not correct")
        phenotype = 'NA'
    return phenotype

def refine_spec_alleles(spec_alleles):
    if len(spec_alleles) == 1:
        spec_alleles = spec_alleles + spec_alleles
    star_alleles  =  []
    print (spec_alleles)
    for i in range(len(spec_alleles)):
        if len(spec_alleles[i].split("*")) > 1:
            spec_alleles[i] = "*" + spec_alleles[i].split("*")[1]
            star_alleles.append("*" + spec_alleles[i].split("*")[1].split('.')[0])
        else:
            star_alleles.append("NA")
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


def output(spec_result, diplotype, detailed_diplotype, phenotype, merge_result_file, tool='spec'):
    f = open(merge_result_file, 'w')
    if tool == 'spec':
        for i in range(1):
            print ('\t'.join(spec_result[i]), file = f)
        print ("#\tCYP2D6 Diplotype:\t"+diplotype, file = f)
        print (f"#\tCYP2D6 Phenotype:\t{phenotype[1]}", file = f)
        print (f"#\tCYP2D6 Activity_score:\t{phenotype[0]}", file = f)
        print ("#\tCYP2D6 Detailed_diplotype\t"+detailed_diplotype, file = f)
        for i in range(1, len(spec_result)):
            print ('\t'.join(spec_result[i]), file = f)
    else:
        print ("#\tCYP2D6 Diplotype:\t"+diplotype, file = f)
        print (f"#\tCYP2D6 Phenotype:\t{phenotype[1]}", file = f)
        print (f"#\tCYP2D6 Activity_score:\t{phenotype[0]}", file = f)
        print ("#\tCYP2D6 Detailed_diplotype:\t"+detailed_diplotype, file = f)
    f.close()

def filter_pangu_hap(pangu_result):
    hap_num_reads = defaultdict(int)

    if len(pangu_result[0]['haplotypes']) > 2:
        over2hap = True
    else:
        over2hap = False

    hap_index = 0
    for hap_dict in pangu_result[0]['haplotypes']:
        hap_reads = 0
        for allele_dict in hap_dict['alleles']:
            allele = allele_dict['call']
            num_reads = allele_dict['num_reads']
            hap_reads += num_reads
        hap_num_reads[hap_index] = hap_reads
        hap_index += 1

    ## sort the alleles by the number of reads
    sorted_alleles = sorted(hap_num_reads.items(), key=lambda x: x[1], reverse=True)
    original_hap = pangu_result[0]['haplotypes'].copy()
    pangu_result[0]['haplotypes'] = []
    pangu_alleles = set()
    diploid_list = []
    for i in range(min(2, len(sorted_alleles))):
        pangu_result[0]['haplotypes'].append(original_hap[sorted_alleles[i][0]])
        diploid_list.append(original_hap[sorted_alleles[i][0]]['call'])
        for j in range(len(original_hap[sorted_alleles[i][0]]['alleles'])):
            pangu_alleles.add(original_hap[sorted_alleles[i][0]]['alleles'][j]['call'])
    diplotype = '/'.join(diploid_list)
    pangu_result[0]['diplotype'] = f"CYP2D6 {diplotype}"
    print ("filtered diplotype", diplotype, pangu_alleles, pangu_result[0]['diplotype'])
    return pangu_result, diplotype, pangu_alleles, over2hap


# refine amplicon output
def refine_amplicon_output(pangu_result, pangu_alleles):
    # print (pangu_result, pangu_alleles)
    allele_num_reads = defaultdict(int)

    for hap_dict in pangu_result[0]['haplotypes']:
        for allele_dict in hap_dict['alleles']:
            allele = allele_dict['call']
            num_reads = allele_dict['num_reads']
            meanCover = allele_dict['meanCover']
            if meanCover == 0:
                continue
            allele_num_reads[allele] += num_reads
    # print (allele_num_reads)
    pangu_alleles = set()
    ## sort the alleles by the number of reads
    sorted_alleles = sorted(allele_num_reads.items(), key=lambda x: x[1], reverse=True)
    if len(sorted_alleles) > 2:
        if sorted_alleles[1][1]/sorted_alleles[0][1] > 0.2:  ## check if homo or hete by reads num
            
            for allele in sorted_alleles[:2]:
                pangu_alleles.add(allele[0])
        else:
            for allele in sorted_alleles[:1]:
                pangu_alleles.add(allele[0])
    else:
        for allele in sorted_alleles:
            pangu_alleles.add(allele[0])
    if len(pangu_alleles) == 1:
        pangu_alleles = list(pangu_alleles) + list(pangu_alleles)
    diplotype = '/'.join(pangu_alleles)
    print (diplotype, pangu_alleles)
    return diplotype, pangu_alleles

def read_star_result(star_result_file):
    ## check if it exists
    if not os.path.exists(star_result_file):
        print (f"WARNING: {star_result_file} does not exist")
        return ['NA', 'NA']
    star_result = []
    with open(star_result_file, 'r') as f:
        f.readline()
        for line in f:
            field = line.strip().split()
            if field[0] == 'cyp2d6':
                star_result.append(field[2])
    print ("star_result", star_result)
    return star_result

if __name__ == "__main__":  

    # prefix = "/mnt/d/HLAPro_backup/Nanopore_optimize/cyp_results/HG00436_1/HG00436_1"
    prefix = sys.argv[1]
    seq_tech = sys.argv[2]
    seq_data = sys.argv[3]
    tool = 'spec' #spec or star

    pangu_result = f"{prefix}_report.json"
    spec_result_file = f"{prefix}.CYP.final.type.result.txt"
    star_result_file = f"{prefix}.Stargazer.type.txt"
    merge_result_file = f"{prefix}.CYP.merge.type.result.txt"
    pangu_result, pangu_alleles = read_pangu_result(pangu_result)
    spec_result, all_spec_result = read_spec_result(spec_result_file)
    if tool != 'spec':
        spec_result = read_star_result(star_result_file)

    if seq_tech == 'amplicon' and len(pangu_result) >= 1:
        diplotype, pangu_alleles = refine_amplicon_output(pangu_result, pangu_alleles)
        detailed_diplotype = 'NA'
    else:
        over2hap = False
        if len(pangu_result) >= 1:
            pangu_result, diplotype, pangu_alleles, over2hap = filter_pangu_hap(pangu_result)
        diplotype, detailed_diplotype = merge_result(pangu_result, pangu_alleles, spec_result, over2hap, tool)
    phenotype = get_phenotype(diplotype)
    output(all_spec_result, diplotype, detailed_diplotype, phenotype, merge_result_file, tool)
