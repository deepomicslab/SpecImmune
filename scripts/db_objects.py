import os

class My_db():

    def __init__(self, args):

        self.gene_class = args["i"]
        self.individual_ref_dir = f"""{args["o"]}/{args["n"]}/individual_ref"""
        if not os.path.exists(self.individual_ref_dir):
            os.makedirs(self.individual_ref_dir)

        if self.gene_class == "HLA":
            self.full_db = f"""{args["db"]}/whole/HLA.full.fasta"""   # 15578 alleles
            self.lite_db = f"""{args["db"]}/whole/HLA.lite.fasta"""   # 6172 alleles
            self.subdir = "whole"

        elif self.gene_class == "KIR":
            # self.lite_db = f"""{args["db"]}/KIR/ref/KIR.extend.select.fasta""" ## 72
            self.lite_db = f"""{args["db"]}/KIR/KIR.full.fasta"""  ## 848
            self.full_db = f"""{args["db"]}/KIR/KIR.lite.fasta"""  ## 848
            self.subdir = "KIR"

        elif self.gene_class == "CYP":
            self.lite_db = f"""{args["db"]}/CYP/ref/CYP.merge.fasta"""  # 1020
            # self.lite_db = f"""{args["db"]}/CYP/ref/CYP.select.fasta"""  # 537
            self.full_db = f"""{args["db"]}/CYP/ref/CYP.merge.fasta"""  # 1020
            self.subdir = "CYP"
        
        elif self.gene_class == "IG_TR":
            self.full_db = f"""{args["db"]}/IG_TR/IG.TR.merge.allele.fasta"""
            self.lite_db = self.full_db
            self.subdir = "IG_TR"
        else:
            print ("wrong gene_class")

        
        self.gene_all_alleles_dir =  f"""{args["db"]}/split_ref/"""  
        self.root = args["db"]

        self.all_alleles = f"""{args["db"]}/whole/merge.fasta"""    ### 17446
        self.allele_len_dict = {}

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

     