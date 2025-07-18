import pandas as pd
from sklearn.preprocessing import StandardScaler
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from collections import defaultdict




def load_matrix():
    df = pd.read_csv("allele_pop_freqs.tsv", sep="\t", header=0)
    ## select top 100 nodes
    # df = df.iloc[:100, :]
    # print (df)
    allele_list = df["Allele"].tolist()
    family_list = df["Family"].tolist()
    ## remove the family column, and covert it to matrix
    df = df.drop(columns=["Family", "Allele"])
    ## remove row name and colunn name

    matrix = df.to_numpy().T
    print (matrix)
    return matrix, allele_list, family_list

def get_gene_precision_matrix():
    X, allele_list, family_list = load_matrix()
    pX_std = StandardScaler().fit_transform(X)

    # from sklearn.covariance import GraphicalLassoCV
    # model = GraphicalLassoCV()
    # model.fit(pX_std)

    from sklearn.covariance import LedoitWolf

    model = LedoitWolf()
    model.fit(pX_std)

    precision = model.precision_


    precision = model.precision_  # Inverse covariance matrix
    print ("Precision matrix shape:", precision.shape)
    print ("Precision matrix:", precision)

    edge_dict = defaultdict(float)
    family_dict = {}
    for i in range(len(allele_list)):
        for j in range(i+1, len(allele_list)):
            gene_1 = allele_list[i].split("*")[0]
            gene_2 = allele_list[j].split("*")[0]
            family_dict[gene_1] = family_list[i]
            family_dict[gene_2] = family_list[j]
            gene_list = sorted([gene_1, gene_2])
            tag = f"{gene_list[0]}_{gene_list[1]}"
            if abs(precision[i, j]) > edge_dict[tag]:  # threshold
                edge_dict[tag] = abs(precision[i, j])
                # edge_dict[tag] = precision[i, j]
    G = nx.Graph()
    for tag, weight in edge_dict.items():
        gene_1, gene_2 = tag.split("_")
        if weight > 0.015:
            G.add_edge(gene_1, gene_2, weight=weight)
            G.add_node(gene_1, label=gene_1, type=family_dict[gene_1])
            G.add_node(gene_2, label=gene_2, type=family_dict[gene_2])


    nx.write_gml(G, "gene_precision_matrix_network.gml")

def get_precision_matrix():
    X, allele_list, family_list = load_matrix()
    pX_std = StandardScaler().fit_transform(X)

    # from sklearn.covariance import GraphicalLassoCV
    # model = GraphicalLassoCV()
    # model.fit(pX_std)

    from sklearn.covariance import LedoitWolf

    model = LedoitWolf()
    model.fit(pX_std)

    precision = model.precision_


    precision = model.precision_  # Inverse covariance matrix
    print ("Precision matrix shape:", precision.shape)
    print ("Precision matrix:", precision)
    threshold = 0.02
    G = nx.Graph()

    node_family_list ={}
    for i in range(len(allele_list)):
        for j in range(i+1, len(allele_list)):
            if abs(precision[i, j]) > threshold:
            # if precision[i, j] < -1 * threshold:
                G.add_edge(allele_list[i], allele_list[j], weight=abs(precision[i, j]))
                node_family_list[allele_list[i]] = family_list[i]
                node_family_list[allele_list[j]] = family_list[j]
    for node in node_family_list:
        G.add_node(node, label=node, type=node_family_list[node])
    # plt.figure(figsize=(12, 12))
    # pos = nx.spring_layout(G, seed=1)  # positions for all nodes
    # nx.draw_networkx_nodes(G, pos, node_size=70, node_color='blue')
    # nx.draw_networkx_edges(G, pos, width=0.5, alpha=0.5)
    # nx.draw_networkx_labels(G, pos, font_size=5, alpha=0.7, font_color="black", font_family="sans-serif")
    # plt.title("Precision Matrix Network")
    # plt.axis("off")
    # plt.savefig("precision_matrix_network.png", format="PNG", dpi=300)
    ## save the graph to gml file
    nx.write_gml(G, "precision_matrix_network.gml")

# matrix = load_matrix()
# print (matrix)
get_precision_matrix()
# get_gene_precision_matrix()
