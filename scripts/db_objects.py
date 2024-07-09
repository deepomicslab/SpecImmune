import os
import re
import numpy as np
import sys



class My_db():

    def __init__(self, args):

        self.gene_class = args["i"]
        
        if "o" in args and "n" in args:
            self.individual_ref_dir = f"""{args["o"]}/{args["n"]}/individual_ref"""
            if not os.path.exists(self.individual_ref_dir):
                os.makedirs(self.individual_ref_dir)

        if self.gene_class == "HLA":
            self.full_db = f"""{args["db"]}/HLA/HLA.full.fasta"""   # 15578 alleles
            self.lite_db = f"""{args["db"]}/HLA/HLA.lite.fasta"""   # 6172 alleles
            self.full_cds_db = f"""{args["db"]}/HLA_CDS/HLA.full.fasta"""
            self.g_group_annotation = f"""{args["db"]}/HLA/hla_nom_g.txt"""
            self.subdir = "HLA"
            # self.subdir = "HLA_CDS" if args["seq_tech"] == "rna" else "HLA"

        elif self.gene_class == "KIR":
            # self.lite_db = f"""{args["db"]}/KIR/ref/KIR.extend.select.fasta""" ## 72
            self.lite_db = f"""{args["db"]}/KIR/KIR.full.fasta"""  ## 848
            self.full_db = f"""{args["db"]}/KIR/KIR.full.fasta"""  ## 848
            self.subdir = "KIR"

        elif self.gene_class == "CYP":
            self.lite_db = f"""{args["db"]}/CYP/CYP.full.fasta"""  
            self.full_db = f"""{args["db"]}/CYP/CYP.full.fasta"""  
            self.subdir = "CYP"
        
        elif self.gene_class == "IG_TR":
            self.full_db = f"""{args["db"]}/IG_TR/IG_TR.fasta"""
            # self.full_db = f"""{args["db"]}/IG_TR/IG.TR.merge.allele.fasta"""
            self.lite_db = self.full_db
            self.subdir = "IG_TR"
            # self.gene_region = f"""{args["db"]}/IG_TR/IG.TR.ref.hap.txt"""
            self.gene_region = f"""{sys.path[0]}/../gene_dist/IG_TR.gene.bed"""
            self.gene_interval_dict, self.order_gene_list = self.get_IG_TR_gene_interval()
            if "hg38" in args:
                self.hg38 = args["hg38"]
                # check if hg38 exists, else exit
                if not os.path.exists(self.hg38):
                    print (f"{self.hg38} does not exist")
                    sys.exit(1)
                # check if hg38 is index by samtools and bwa, else index it
                if not os.path.exists(f"{self.hg38}.fai"):
                    os.system(f"samtools faidx {self.hg38}")
                if not os.path.exists(f"{self.hg38}.bwt"):
                    os.system(f"bwa index {self.hg38}")
            else:
                print ("WARNING: hg38 not provided, cannot type IG_TR genes")
                # sys.exit(1)
        else:
            print ("wrong gene_class")

        
        self.gene_all_alleles_dir =  f"""{args["db"]}/split_ref/"""  
        self.root = args["db"]

        self.all_alleles = f"""{args["db"]}/whole/merge.fasta"""    ### 17446
        self.allele_len_dict = {}
        self.version_info = "# version:  N/A"
        self.get_db_version(args)

        self.gene_min_len = {}
        self.cal_gene_len()

    def get_db_version(self, args):
        g_file = "%s/%s/release_version.txt"%(args['db'], self.subdir)

        ## if the release_version.txt does not exist, return "N/A"
        if not os.path.exists(g_file):
            print (f"release_version.txt does not exist in {g_file}")
            pass
        else:
            for line in open(g_file):
                if re.search("# version:", line):
                    self.version_info = line.strip()

    def get_gene_alleles(self, gene):
        ### record representative allele of each gene
        # return self.gene_all_alleles_dir + f"/{gene}.fasta"
        return self.individual_ref_dir + f"/{gene}/{gene}.fasta"

    def get_gene_alleles_2ref(self, gene, ref_idx):
        ### record representative allele of each gene
        return self.individual_ref_dir + f"/{gene}/{gene}.{ref_idx+1}.fasta"
    
    def get_gene_alleles_ref(self, gene):
        ### record representative allele of each gene
        return self.individual_ref_dir + f"/{gene}/{gene}.fasta"

    def get_gene_all_alleles(self, gene):
        ### record all alleles of each gene

        # if self.gene_class == "HLA" or self.gene_class == "CYP":
        #     ref = f"""{self.root}/{self.gene_class}/ref/split/{gene}.fasta"""
        # elif self.gene_class == "KIR":
        #     ref = f"""{self.root}/KIR/ref/split/{gene}.exon.fasta"""
        # else:
        #     print ("wrong gene_class")
        # ref = f"""{self.root}/whole/{gene}.fasta"""
        ref = f"""{self.root}/{self.subdir}/{gene}/{gene}.fasta"""
        return ref

    def get_allele_length(self, gene):
        fai_file = f"""{self.root}/{self.subdir}/{gene}/{gene}.fasta.fai"""
        # fai_file = f"""{self.root}/whole/{gene}/{gene}.fasta.fai"""
        ## check if faidx file exists
        if not os.path.exists(fai_file):
            print (fai_file)
            print ("fai file does not exist")
            return
        f = open(fai_file, "r")
        for line in f:
            length = int(line.split("\t")[1])
            allele = line.split("\t")[0]
            self.allele_len_dict[allele] = length
        # print (self.allele_len_dict)
        f.close()

    def get_blast_index(self, gene):

        ref = f"""{self.root}/{self.subdir}/{gene}/{gene}"""


        return ref

    def cal_gene_len(self):
        gene_length_dict = {}
        db_dir = f"{self.root}/{self.subdir}"
        # for each dir in the db_dir
        for gene in os.listdir(db_dir):
            gene_dir = os.path.join(db_dir, gene)
            if os.path.isdir(gene_dir):
                for file in os.listdir(gene_dir):
                    if file.endswith(".fasta"):
                        fai_file = f"{self.full_db}.fai"
                        if not os.path.exists(fai_file):
                            # print (fai_file)
                            print (fai_file, "fai file does not exist")
                            return
                        gene_length_dict[gene] = []
                        f = open(fai_file, "r")
                        for line in f:
                            length = int(line.split("\t")[1])
                            gene_length_dict[gene].append(length)
                        f.close()
        for gene in gene_length_dict:
            # self.gene_min_len[gene] = min(gene_length_dict[gene])
            self.gene_min_len[gene] = round(np.mean(gene_length_dict[gene]))
        # print (self.gene_min_len["MICA"])

    def get_IG_TR_gene_interval(self):
        gene_list = []
        interval_dict = {}
        with open(self.gene_region, "r") as f:
            for line in f:
                line = line.strip().split()
                gene_list.append(line[0])
                if line[1] not in interval_dict:
                    interval_dict[line[1]] = {}
                interval_dict[line[1]][line[0]] = (int(line[2]), int(line[3]))
        # print ("intervali_dict",interval_dict)
        order_gene_list = []
        for genome in interval_dict:
            # sort the gene on it by their start position
            interval_dict[genome] = dict(sorted(interval_dict[genome].items(), key=lambda item: item[1][0]))
            for gene in interval_dict[genome]:
                order_gene_list.append(gene)

        return interval_dict, order_gene_list


"""
def read_G_annotation(db):
    g_file = f"{db}/HLA/hla_nom_g.txt"
    G_annotation_dict = {}
    i = 0
    for line in open(g_file):
        if re.search("# version:", line):
            version_info = line
        if line[0] == "#":
            continue
        # line.replace(":","_", 10000)
        line = re.sub(":","_",line)
        array = line.strip().split(";")
        gene = array[0][:-1]
        
        if len(array[-1]) == 0:
            
            g_name = gene + "_" + array[-2]
            # print (g_name)
            G_annotation_dict[g_name] = g_name
        else:
            g_name = gene + "_" + array[-1]
            alleles = array[-2].split("/")
            for each in alleles:
                each = gene + "_" + each
                G_annotation_dict[each] = g_name
        # print (array, g_name)
        # print (G_annotation_dict)
        # print (len(array))
        # if i > 2:
        #     break
        i += 1
    return G_annotation_dict, version_info
"""