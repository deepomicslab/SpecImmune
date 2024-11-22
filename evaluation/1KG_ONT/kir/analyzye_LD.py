import os
import pandas as pd

def count_line_number_of_file(file):
    with open(file, 'r') as f:
        i = 0
        for line in f:
            i += 1
    return i 

def read_LD_values(raw_dir, indir,outfile, sample_num_threshold=500):
    ## for each file in the indir
    ## read the LD values and store them in a df
    ## return the df
    data = []
    gene_set = set()
    for file in os.listdir(indir):
        if file.endswith(".csv"):
            # print(file)
            with open(indir + file, 'r') as f:
                i = 0
                for line in f:
                    if i == 1:
                        D = line.strip().replace('"', '')
                    if i == 2:
                        Wn = line.strip().replace('"', '')
                    if i == 3:
                        Wab = line.strip().replace('"', '')
                    if i == 4:
                        Wba = line.strip().replace('"', '')
                    if i == 5:
                        hap_num = line.strip().replace('"', '')
                    i += 1
            sample_num = count_line_number_of_file(raw_dir + file) -1
            if sample_num < sample_num_threshold:
                continue
            if Wn == "Inf" or Wn == "NaN":
                continue
            field = file.split("_") 
            gene1 = field[0]
            gene2 = field[1]
            min_w = min(float(Wab), float(Wba))
            data.append([gene1, gene2, D, Wn, hap_num, sample_num, Wab, min_w])
            data.append([gene2, gene1, D, Wn, hap_num, sample_num, Wba, min_w])
            gene_set.add(gene1)
            gene_set.add(gene2)

            # print (gene1, gene2, D, Wn)
            # break
    for gene in gene_set:
        data.append([gene, gene, 1, 1, 1000, 1000, 1, 1])
    df = pd.DataFrame(data, columns = ['gene1', 'gene2', 'D', 'Wn', 'hap_num','sample_num','Wab', 'min_w'])
    df.to_csv(outfile, index = False)


def graph(outfile):
    # df3 <- df[df$min_w > 0.2,]
    # ## covert df3 to a matrix

    # df_matrix3 <- reshape2::acast(df3, gene1 ~ gene2, value.var = "min_w", fill = "min_w")
    # ### save the matrix to csv file
    # write.csv(df_matrix3, "minALD_matrix.tsv")
    ## use network x read the graph
    import networkx as nx
    import numpy as np
    import pandas as pd

    ## read the csv
    df = pd.read_csv(outfile)
    df3 = df[df['min_w'] > 0]
    # df3 = df[df['min_w'] > 0.2]
    ## fill the matrix with zero is NA
    df_matrix3 = df3.pivot(index='gene1', columns='gene2', values='min_w')
    ## convert NaN to 0
    df_matrix3 = df_matrix3.fillna(0)
    print (df_matrix3)
    ## output the matrix to a tsv file
    df_matrix3.to_csv("tree/minALD_matrix.tsv", sep = "\t")
    ## load the matrix to a graph using networkx
    G = nx.from_pandas_adjacency(df_matrix3)
    ## compute the importance of centrality of each node, and sort the node by centrality, and output to a csv file
    centrality = nx.eigenvector_centrality(G)
    centrality = sorted(centrality.items(), key=lambda x:x[1], reverse=True)
    with open("tree/centrality.tsv", 'w') as f:
        for node in centrality:
            print (node[0], node[1], file = f)
    ## compute page rank score of each node, and sort the node by score, and output to a csv file
    page_rank = nx.pagerank(G)
    page_rank = sorted(page_rank.items(), key=lambda x:x[1], reverse=True)
    with open("tree/page_rank.tsv", 'w') as f:
        for node in page_rank:
            print (node[0], node[1], file = f)
    ### compute the betweeness of each node, and sort the node by betweeness, and print the top 10 nodes    
    betweeness = nx.betweenness_centrality(G)
    betweeness = sorted(betweeness.items(), key=lambda x:x[1], reverse=True)
    # with open("tree/betweeness.tsv", 'w') as f:
    for node in betweeness[:5]:
        print ("betweeness", node[0], node[1])
    ## compute the clustering coefficient of each node, and sort the node by clustering coefficient, and print the top 10 nodes
    clustering = nx.clustering(G)
    clustering = sorted(clustering.items(), key=lambda x:x[1], reverse=True)
    # with open("tree/clustering.tsv", 'w') as f:
    for node in clustering[:5]:
        print ("clustering coefficient", node[0], node[1])
    ## compute the degree of each node, and sort the node by degree, and print the top 10 nodes
    degree = nx.degree_centrality(G)
    degree = sorted(degree.items(), key=lambda x:x[1], reverse=True)
    # with open("tree/degree.tsv", 'w') as f:
    for node in degree[:5]:
        print ("degree", node[0], node[1])



# indir = "/mnt/d/HLAPro_backup/Nanopore_optimize/1kgp_analysis/kir_LD_result/"
# outfile = "kir_LD_values.csv"
# read_LD_values(indir,outfile)

# indir = "/mnt/d/HLAPro_backup/Nanopore_optimize/1kgp_analysis/hla_LD_result/"
# outfile = "hla_LD_values.csv"
# read_LD_values(indir,outfile)

# raw_dir = "/mnt/d/HLAPro_backup/Nanopore_optimize/1kgp_analysis/hla_kir_LD/"
# indir = "/mnt/d/HLAPro_backup/Nanopore_optimize/1kgp_analysis/hla_kir_LD_result/"
# outfile = "hla_kir_LD_values.csv"
# read_LD_values(raw_dir, indir,outfile)

# raw_dir = "/mnt/d/HLAPro_backup/Nanopore_optimize/1kgp_analysis/hla_kir_cyp_LD/"
# indir = "/mnt/d/HLAPro_backup/Nanopore_optimize/1kgp_analysis/hla_kir_cyp_LD_result/"
# outfile = "hla_kir_cyp_LD_values.csv"
# read_LD_values(raw_dir, indir,outfile)

# print ("#####################")
# raw_dir = "/mnt/d/HLAPro_backup/Nanopore_optimize/1kgp_analysis/hla_kir_cyp_vdj_LD/"
# indir = "/mnt/d/HLAPro_backup/Nanopore_optimize/1kgp_analysis/hla_kir_cyp_vdj_LD_result/"
# outfile = "hla_kir_cyp_vdj_LD_values.csv"
# read_LD_values(raw_dir, indir,outfile)

# print ("#####################")
raw_dir = "/mnt/d/HLAPro_backup/Nanopore_optimize/1kgp_analysis/hla_kir_cyp_vdj_LD2/"
indir = "/mnt/d/HLAPro_backup/Nanopore_optimize/1kgp_analysis/hla_kir_cyp_vdj_LD_result2/"
outfile = "hla_kir_cyp_vdj_LD_values2.csv"
# read_LD_values(raw_dir, indir,outfile)

graph(outfile)
