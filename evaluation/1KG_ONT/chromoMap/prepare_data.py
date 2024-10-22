import pandas as pd

## load the bed file
def read_bed(bed_file, gene_class='HLA'):
    data = []
    f = open(bed_file, 'r')
    for line in f:
        if line.startswith('#'):
            continue
        else:
            field = line.strip().split()
            if field[1] == '-':
                continue
            # if field[1] != 'chr14':
            #     continue
            if field[0] == 'HFE':
                continue
            chr_id = field[1]
            start = int(field[2])
            end = field[2]
            if chr_id == 'chr14':
                if start < 100363198:
                    chr_id = 'chr14_part1'
                else:
                    chr_id = 'chr14_part2'
            if chr_id == 'chr7':
                if start < 100799177:
                    chr_id = 'chr7_part1'
                else:
                    chr_id = 'chr7_part2'
            if chr_id == 'chr22':
                if start < 25000000:
                    chr_id = 'chr22_part1'
                else:
                    chr_id = 'chr22_part2'
            
            if gene_class == "CYP" and field[0] != "CYP2D6":
                continue
            field[1] = chr_id

            data.append(field)
    return data

def handle_data(data):
    chr_data = {}
    anno_data = {}
    for line in data:
        # print (line)
        gene, chr_id, start, end, anno = line
            
        start = int(start)
        end = int(end)

        if chr_id not in chr_data:
            chr_data[chr_id] = [min(start, end), max(start, end)]
        if min(start, end) < chr_data[chr_id][0]:
            chr_data[chr_id][0] = min(start, end)
        if max(start, end) > chr_data[chr_id][1]:
            chr_data[chr_id][1] = max(start, end)
        anno_data[gene] = line
    # output(chr_data, anno_data)
    personal_output(chr_data, anno_data, sample_info)

def output(chr_data, anno_data):
    
    for chr_id in chr_data:
        print (chr_id, chr_data[chr_id])
        chr_file = f"/mnt/d/R_script_files/mengyao/test_chr_{chr_id}.txt"
        chr_out = open(chr_file, 'w')
        print(chr_id, chr_data[chr_id][0], chr_data[chr_id][1], sep = '\t', file=chr_out)
        chr_out.close()

        anno_file = f"/mnt/d/R_script_files/mengyao/test_anno_{chr_id}.txt"
        anno_out = open(anno_file, 'w')
        for gene in anno_data:
            if anno_data[gene][1] == chr_id:
                print("\t".join(anno_data[gene][:4]+['1']), file=anno_out)
        anno_out.close()

def personal_output(chr_data, anno_data, sample_info):
    
    for chr_id in chr_data:
        print (chr_id, chr_data[chr_id])
        chr_file = f"/mnt/d/R_script_files/mengyao/test_chr_{chr_id}.txt"
        chr_out = open(chr_file, 'w')
        print(chr_id, chr_data[chr_id][0], chr_data[chr_id][1], sep = '\t', file=chr_out)
        chr_out.close()

        anno_file = f"/mnt/d/R_script_files/mengyao/test_anno_{chr_id}.txt"
        anno_out = open(anno_file, 'w')
        for gene in anno_data:
            if gene not in sample_info:
                continue
            else:
                genotype = sample_info[gene]
            if anno_data[gene][1] == chr_id:
                anno_data[gene][0] += "*"+genotype
                print("\t".join(anno_data[gene][:4]+['1']), file=anno_out)
        anno_out.close()


def read_all():
    all_data = []
    all_data += read_bed("../../../gene_dist/HLA.gene.bed")
    all_data += read_bed("../../../gene_dist/KIR.gene.bed")
    all_data += read_bed("../../../gene_dist/CYP.gene.bed", 'CYP')
    all_data += read_bed("../../../gene_dist/IG_TR.gene.bed")
    handle_data(all_data)

def get_type():
    typed_loci_num = {}
    f = open(all_types, 'r')
    for line in f:
        field = line.strip().split(',')
        sample = field[0]
        sample_info = {}
        for i in range(1, len(field)):
            gene_info = field[i].split("*")
            gene = gene_info[0]
            genotype = gene_info[1]
            # print (sample, chr_file, anno_file)
            if gene not in sample_info:
                sample_info[gene] = genotype
            else:
                sample_info[gene] += "/"+genotype
        # print (sample, len(sample_info))
        typed_loci_num[sample] = len(sample_info)
        if sample == 'HG00377':
            # print (sample_info['CYP2D6'])
            return sample_info
    ## sort the samples by the number of typed loci
    sorted_samples = sorted(typed_loci_num.items(), key=lambda x:x[1], reverse=True)
    print (sorted_samples[:5])



chr_file = "/mnt/d/R_script_files/mengyao/test_chr.txt"
anno_file = "/mnt/d/R_script_files/mengyao/test_anno.txt"

# hla = '/home/wangshuai/softwares/SpecLong/evaluation/1KG_ONT/hla/speclong_res_merged_samples.csv'
# kir = '/home/wangshuai/softwares/SpecLong/evaluation/1KG_ONT/kir/merged_samples.csv'
# cyp = '/home/wangshuai/softwares/SpecLong/evaluation/1KG_ONT/cyp/cyp_1k_all.csv'
# vdj = '/home/wangshuai/softwares/SpecLong/evaluation/1KG_ONT/ig_tr/merged_samples.ig_tr.csv'
all_types = "sample_all_loci_type.csv"
sample_info = get_type()

read_all()