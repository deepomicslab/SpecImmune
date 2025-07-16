import pandas as pd
import networkx as nx
# Draw the graph
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from networkx.algorithms.community import label_propagation_communities


def profile_graph(G, family_dict):
    ## report the top 10 nodes with the highest degree
    top_nodes = sorted(G.degree, key=lambda x: x[1], reverse=True)[:5]
    print("Top 5 nodes by degree:")
    for node, degree in top_nodes:
        print(f"{node}: {degree}")

    ## report the top 5 nodes with the highest betweenness centrality
    top_betweenness = sorted(nx.betweenness_centrality(G).items(), key=lambda x: x[1], reverse=True)[:5]
    print("\nTop 5 nodes by betweenness centrality:")
    for node, centrality in top_betweenness:
        print(f"{node}: {centrality:.4f}")

    ## report the top 5 nodes with the highest page rank
    top_page_rank = sorted(nx.pagerank(G).items(), key=lambda x: x[1], reverse=True)[:5]
    print("\nTop 5 nodes by PageRank:")
    for node, rank in top_page_rank:
        print(f"{node}: {rank:.4f}")

    ## report the density of the graph and node number
    density = nx.density(G)
    print(f"\nGraph density: {density:.4f}")
    print(f"Number of nodes: {G.number_of_nodes()}")
    print(f"Number of edges: {G.number_of_edges()}")
    ## count how many nodes are in each family
    family_count = {}
    for node in G.nodes():
        family = family_dict.get(node, "Unknown")
        if family not in family_count:
            family_count[family] = 0
        family_count[family] += 1

    print("\nNumber of nodes in each family:")
    for family, count in family_count.items():
        print(f"{family}: {count}")
    
    ## if the graph is connected
    if nx.is_connected(G):
        print("\nThe graph is connected.")
    else:
        print("\nThe graph is not connected.")


def cluster(G, family_dict):

    # import community as community_louvain
    # partition = community_louvain.best_partition(G)
    # print("\nCommunity structure (Louvain method):")
    # community_count = {}
    # for node, community in partition.items():
    #     if community not in community_count:
    #         community_count[community] = 0
    #     community_count[community] += 1

    # for community, count in community_count.items():
    #     print(f"Community {community}: {count} nodes")
    
    # ## print the node of each community
    # print("\nNodes in each community:")
    # community_nodes = {}
    # for node, community in partition.items():
    #     if community not in community_nodes:
    #         community_nodes[community] = []
    #     community_nodes[community].append(node) 
    # for community, nodes in community_nodes.items():
    #     print(f"Community {community}: {', '.join(sorted(nodes))}")

    
    communities = label_propagation_communities(G)
    print("\nCommunity structure (Label Propagation method):")
    for i, community in enumerate(communities):
        print(f"Community {i}: {len(community)} nodes")
    ## print the nodes in each community
    print("\nNodes in each community (Label Propagation):")
    for i, community in enumerate(communities):
        community_nodes = sorted(community)
        print(f"Community {i}: {', '.join(community_nodes)}")
    ### save the community structure to a CSV file
    data = []
    for i, community in enumerate(communities):
        if i == 0:
            first = community
        elif i == 1:
            second = community
        for node in community:
            data.append({"Community": i, "Node": node, "Family": family_dict.get(node, "Unknown")})
    partition_df = pd.DataFrame(data)
    partition_df.to_csv("community_structure.csv", index=False)

    data_list = []
    for u, v, data in G.edges(data=True):
        if u in first and v in second:  # Example condition to filter edges
            data_list.append([u, v, data['weight']])
        elif v in first and u in second:
            data_list.append([v, u, data['weight']])
    edge_df = pd.DataFrame(data_list, columns=["Community_0", "Community_1", "Weight"])
    edge_df.to_csv("community_edges.csv", index=False)


def plot(G, family_dict):
    ## group the nodes by family
    family_groups = {}
    for node in G.nodes():
        family = family_dict.get(node, "Unknown")
        if family not in family_groups:
            family_groups[family] = []
        family_groups[family].append(node)



    plt.figure(figsize=(12, 12))
    pos = nx.spring_layout(G, seed=1)  # positions for all nodes

    # Assign a color to each family
    families = list(family_groups.keys())
    colors = list(mcolors.TABLEAU_COLORS.values())
    family_color_map = {fam: colors[i % len(colors)] for i, fam in enumerate(families)}

    # Draw nodes with colors by family
    for fam, nodes in family_groups.items():
        nx.draw_networkx_nodes(G, pos, nodelist=nodes, node_size=70, node_color=family_color_map[fam], label=fam)

    nx.draw_networkx_edges(G, pos, width=0.5, alpha=0.5)
    nx.draw_networkx_labels(G, pos, font_size=5, alpha=0.7, font_color="black", font_family="sans-serif")
    plt.legend(title="Family")
    ## save the plot
    plt.title("Co-evolution Network")
    plt.axis("off")
    plt.savefig("HLA_KIR_co_evolution_network.png", format="PNG", dpi=300)
    plt.show()
    plt.close()

def load_graph():
    ## define a graph
    G = nx.Graph()
    corr_csv = "co_evo_pearson_results.csv"
    df = pd.read_csv(corr_csv)
    IG_TCR = ["TCR", "IG"]
    # df = df[df["corr_tag"] == "HLA_vs_KIR"]
    pair_dict = {}
    family_dict = {}
    for i, row in df.iterrows():
        gene_1 = row["HLA_allele"].split("*")[0]
        # print (i, row["KIR_allele"])
        gene_2 = row["KIR_allele"].split("*")[0]
        r = row["Pearson_r"]
        compare_tag = row["corr_tag"]
        family_list = compare_tag.split("_vs_")
        family_1 = family_list[0]
        family_2 = family_list[1]

        family_dict[gene_1] = family_1
        family_dict[gene_2] = family_2
        G.add_node(gene_1, label=family_dict[gene_1], type=family_dict[gene_1])
        G.add_node(gene_2, label=family_dict[gene_2], type=family_dict[gene_2])

        tag = f"{gene_1}_{gene_2}"
        if tag not in pair_dict:
            pair_dict[tag] = []
        pair_dict[tag].append(abs(r))

    for tag, r_list in pair_dict.items():
        # print (tag, max(r_list), r_list)
        gene_1, gene_2 = tag.split("_")
        G.add_edge(gene_1, gene_2, weight=max(r_list))
    nx.write_gml(G, "co_evolution_network.gml")
    return G, family_dict

if __name__ == "__main__":
    G, family_dict = load_graph()
    # profile_graph(G, family_dict)
    # plot(G, family_dict)
    cluster(G, family_dict)


