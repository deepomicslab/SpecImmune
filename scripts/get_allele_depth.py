from collections import defaultdict
import numpy as np

class Get_depth():

    def __init__(self, depth_file):
        self.depth_file = depth_file
        self.depth_dict = {}

    def record_depth(self):
        f = open(self.depth_file)
        for line in f:
            array = line.strip().split()

            allele = array[0]

            gene = allele.split("*")[0]
            
            depth = int(array[2])
            # print (gene, allele, depth)
            if gene not in self.depth_dict:
                self.depth_dict[gene] = {}
            if allele not in self.depth_dict[gene]:
                self.depth_dict[gene][allele] = []

            self.depth_dict[gene][allele].append(depth)
    
    def select(self, sort_depth_file, max_allele_num=10):
        f = open(sort_depth_file, 'w')
        print ("Gene\tAllele\tDepth\tAllele_length", file = f)

        record_candidate_alleles = defaultdict(set)
        record_allele_length = {}
        for gene in self.depth_dict:
            # if gene not in gene_list:
            #     continue
            record_allele_depth = {}
            
            record_allele_info = {}
            for allele in self.depth_dict[gene]:
                mean_depth = np.mean(self.depth_dict[gene][allele])
                median_depth = np.median(self.depth_dict[gene][allele])

                over_0 = 0
                for e in self.depth_dict[gene][allele]:
                    if e > 0:
                       over_0 += 1
                coverage =  float(over_0)/len(self.depth_dict[gene][allele])
                record_allele_length[allele] = len(self.depth_dict[gene][allele])

                record_allele_depth[allele] = mean_depth
                record_allele_info[allele] = [mean_depth, median_depth, coverage]

                # print (allele, mean_depth, median_depth, coverage)
            sorted_dict = sorted(record_allele_depth.items(), key=lambda x: x[1], reverse=True)
            for i in range(min([max_allele_num, len(sorted_dict)])):
                print (sorted_dict[i][0], round(sorted_dict[i][1],2), record_allele_length[sorted_dict[i][0]] , sep = "\t")
                record_candidate_alleles[gene].add(sorted_dict[i][0])

            ## output
            for i in range(len(sorted_dict)):
                print (gene, sorted_dict[i][0], round(sorted_dict[i][1],2), record_allele_length[sorted_dict[i][0]], sep = "\t", file = f)
        f.close()

        return record_candidate_alleles, record_allele_length