import sys


class Bed_db():

    def __init__(self):

        self.gene_file =  f"{sys.path[0]}/../gene_dist//IG_TR.gene.bed"   #sys.argv[1]
        self.segment_bed = f"{sys.path[0]}/../gene_dist/IG_TR.segment.bed"
        self.lite_gene_file =  f"{sys.path[0]}/../gene_dist//IG_TR.gene.lite.bed" 
    
    def get_gene_interval(self, file):
        interval_dict = {}
        with open(file, "r") as f:
            for line in f:
                line = line.strip().split()
                gene = line[0]
                chrom = line[1]
                start = int(line[2])
                end = int(line[3])

                interval_dict[gene] = f"{chrom}:{start}-{end}"
        return interval_dict
    
    def get_hg38_gene_interval(self):
        hg38_gene_info = {}
        with open(self.gene_file, "r") as f:
            for line in f:
                line = line.strip().split()
                gene = line[0]
                chrom = line[1]
                start = int(line[2])
                end = int(line[3])

                hg38_gene_info[gene] = [chrom, start, end, end - start]
        return hg38_gene_info