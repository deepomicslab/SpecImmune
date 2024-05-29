

class My_db():

    def __init__(self, args):
        self.gene_class = args["i"]
        if self.gene_class == "HLA":
            self.lite_db = f"""{args["db"]}/HLA/ref/HLA.extend.fasta"""
        elif self.gene_class == "KIR":
            self.lite_db = f"""{args["db"]}/KIR/ref/KIR.extend.select.fasta"""
            #self.db = f"{sys.path[0]}/../db/KIR/ref/KIR.extend.fasta"
        elif self.gene_class == "CYP":
            self.lite_db = f"""{args["db"]}/CYP/ref/CYP.merge.fasta"""
        else:
            print ("wrong gene_class")

        
        self.gene_all_alleles_dir =  f"""{args["db"]}/split_ref/"""  
        self.root = args["db"]

    def get_gene_alleles(self, gene):
        ### record representative allele of each gene
        return self.gene_all_alleles_dir + f"/{gene}.fasta"


    def get_gene_all_alleles(self, gene):
        ### record all alleles of each gene

        if self.gene_class == "HLA" or self.gene_class == "CYP":
            ref = f"""{self.root}/{self.gene_class}/ref/split/{gene}.fasta"""
        elif self.gene_class == "KIR":
            ref = f"""{self.root}/KIR/ref/split/{gene}.exon.fasta"""
        else:
            print ("wrong gene_class")
        
        return ref
     