import os

class My_db():

    def __init__(self, args):


        # self.db = f"{sys.path[0]}/../db/ref/hla_gen.format.filter.extend.DRB.no26789.fasta"
        # self.db = f"{sys.path[0]}/../db/ref/hla_gen.format.filter.extend.DRB.no26789.v2.fasta"
        # self.db = f"""{args["db"]}/ref/hla_gen.format.filter.extend.DRB.no26789.fasta"""

        self.gene_class = args["i"]
        self.individual_ref_dir = f"""{args["o"]}/{args["n"]}/individual_ref"""
        if not os.path.exists(self.individual_ref_dir):
            os.makedirs(self.individual_ref_dir)

        if self.gene_class == "HLA":
            self.full_db = f"""{args["db"]}/whole/HLA.full.fasta"""   # 15578 alleles
            self.lite_db = f"""{args["db"]}/whole/HLA.lite.fasta"""   # 6172 alleles
            # self.lite_db = f"/mnt/d/HLAPro_backup/Nanopore_optimize/SpecHLA/db/ref/hla_gen.format.filter.extend.DRB.no26789.fasta"

        elif self.gene_class == "KIR":
            # self.lite_db = f"""{args["db"]}/KIR/ref/KIR.extend.select.fasta""" ## 72
            self.lite_db = f"""{args["db"]}/KIR/ref/KIR.extend.fasta"""  ## 848
            self.full_db = f"""{args["db"]}/KIR/ref/KIR.extend.fasta"""  ## 848

        elif self.gene_class == "CYP":
            self.lite_db = f"""{args["db"]}/CYP/ref/CYP.merge.fasta"""  # 1020
            # self.lite_db = f"""{args["db"]}/CYP/ref/CYP.select.fasta"""  # 537
            self.full_db = f"""{args["db"]}/CYP/ref/CYP.merge.fasta"""  # 1020
        
        elif self.gene_class == "IG_TR":
            self.full_db = f"""{args["db"]}/IG_TR/IG.TR.merge.allele.fasta"""
            self.lite_db = self.full_db

        else:
            print ("wrong gene_class")

        
        self.gene_all_alleles_dir =  f"""{args["db"]}/split_ref/"""  
        self.root = args["db"]

        self.all_alleles = f"""{args["db"]}/whole/merge.fasta"""    ### 17446

    def get_gene_alleles(self, gene):
        ### record representative allele of each gene
        # return self.gene_all_alleles_dir + f"/{gene}.fasta"
        return self.individual_ref_dir + f"/{gene}/{gene}.fasta"
    
    def get_gene_alleles_2ref(self, gene, ref_idx):
        ### record representative allele of each gene
        return self.individual_ref_dir + f"/{gene}/{gene}.{ref_idx+1}.fasta"


    def get_gene_all_alleles(self, gene):
        ### record all alleles of each gene

        # if self.gene_class == "HLA" or self.gene_class == "CYP":
        #     ref = f"""{self.root}/{self.gene_class}/ref/split/{gene}.fasta"""
        # elif self.gene_class == "KIR":
        #     ref = f"""{self.root}/KIR/ref/split/{gene}.exon.fasta"""
        # else:
        #     print ("wrong gene_class")
        # ref = f"""{self.root}/whole/{gene}.fasta"""
        ref = f"""{self.root}/whole/{gene}/{gene}.fasta"""
        # ref = f"""{self.root}/clean_whole/{gene}.fasta"""
        return ref

     