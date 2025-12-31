import pandas as pd
from sklearn.preprocessing import StandardScaler
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from collections import defaultdict
from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression
# from sklearn.covariance import GraphicalLassoCV
# model = GraphicalLassoCV()
# model.fit(pX_std)

from sklearn.covariance import LedoitWolf




def load_matrix():
    """Load raw allele frequencies and apply PC adjustment"""
    df = pd.read_csv("allele_pop_freqs.tsv", sep="\t", header=0)
    ## select top 100 nodes
    # df = df.iloc[:100, :]
    # print (df)
    allele_list = df["Allele"].tolist()
    family_list = df["Family"].tolist()
    ## remove the family column, and covert it to matrix
    df = df.drop(columns=["Family", "Allele"])
    ## remove row name and column name

    matrix = df.to_numpy().T  # Shape: (n_populations, n_alleles)
    
    # Compute PCs from all allele frequencies (following co_evo_pearson.py approach)
    print("Computing PCs for population structure adjustment...")
    pca = PCA(n_components=min(5, matrix.shape[0], matrix.shape[1]))
    pcs = pca.fit_transform(matrix)
    print(f"PCs computed: {pcs.shape}, variance explained: {pca.explained_variance_ratio_}")
    
    # Regress out PCs from each allele frequency
    print("Regressing out PCs from allele frequencies...")
    matrix_adjusted = np.zeros_like(matrix)
    for i in range(matrix.shape[1]):  # For each allele
        freq = matrix[:, i].reshape(-1, 1)
        reg = LinearRegression().fit(pcs, freq)
        residuals = freq - reg.predict(pcs)
        matrix_adjusted[:, i] = residuals.flatten()
    
    print("PC-adjusted allele frequencies computed")
    print(f"Original matrix shape: {matrix.shape}")
    print(f"Adjusted matrix shape: {matrix_adjusted.shape}")
    
    return matrix_adjusted, allele_list, family_list

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

def estimate_cutoff(precision, allele_list):
    ## count the number of edges in the precision matrix using different thresholds
    data = []
    for threshold in np.arange(0.01, 0.05, 0.005):
        edge_count = 0
        for i in range(len(allele_list)):
            for j in range(i+1, len(allele_list)):
                if abs(precision[i, j]) > threshold:
                    edge_count += 1
        data.append((threshold, edge_count))
    df = pd.DataFrame(data, columns=["Threshold", "Edge Count"])
    ## plot the data
    plt.figure(figsize=(8, 5))
    plt.plot(df["Threshold"], df["Edge Count"], marker='o')
    plt.axvline(x=0.02, color='r', linestyle='--', label='Threshold = 0.02', alpha=0.7)
    plt.xlabel("Threshold")
    plt.ylabel("Edge Count")
    plt.title("Sparsity Path: Network Size vs Threshold")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig("plots/edge_count_vs_threshold.pdf")
    plt.close()
    
    ## plot the distribution of the precision matrix, and plot a line at 0.02
    plt.figure(figsize=(8, 5))
    plt.hist(precision.flatten(), bins=100, density=True)
    plt.axvline(x=0.02, color='r', linestyle='--', label='Threshold = 0.02')
    plt.xlabel("Precision Value")
    plt.ylabel("Density")
    plt.title("Distribution of Precision Matrix Values")
    plt.legend()
    plt.savefig("plots/precision_matrix_distribution.pdf")
    plt.close()
    
    ## count the how proportion of values are above the threshold
    above_threshold = np.sum(np.abs(precision) > 0.02)
    total_values = precision.size
    proportion_above_threshold = above_threshold / total_values
    print(f"Proportion of values above threshold 0.02: {proportion_above_threshold:.4f}")
    
    # Save sparsity path data
    df.to_csv("plots/sparsity_path.csv", index=False)
    return df


def stability_selection(X, allele_list, family_list, n_bootstrap=100, subsample_ratio=0.8, threshold=0.02):
    """
    Perform stability selection through subsampling to assess edge robustness.
    
    Parameters
    ----------
    X : ndarray
        PC-adjusted allele frequency matrix (n_populations x n_alleles)
    allele_list : list
        List of allele names
    family_list : list
        List of family names for each allele
    n_bootstrap : int
        Number of bootstrap/subsampling iterations
    subsample_ratio : float
        Proportion of populations to sample in each iteration
    threshold : float
        Edge weight threshold for inclusion
    
    Returns
    -------
    stability_scores : dict
        Dictionary mapping edge (i,j) to its stability score (proportion of times selected)
    """
    print(f"\nPerforming stability selection with {n_bootstrap} bootstrap iterations...")
    print(f"Subsample ratio: {subsample_ratio}, Threshold: {threshold}")
    
    n_populations = X.shape[0]
    n_alleles = X.shape[1]
    n_subsample = int(n_populations * subsample_ratio)
    
    # Track edge selections across bootstrap samples
    edge_selections = defaultdict(int)
    
    for boot_iter in range(n_bootstrap):
        if (boot_iter + 1) % 20 == 0:
            print(f"  Bootstrap iteration {boot_iter + 1}/{n_bootstrap}")
        
        # Subsample populations with replacement
        subsample_idx = np.random.choice(n_populations, size=n_subsample, replace=True)
        X_sub = X[subsample_idx, :]
        
        # Standardize and fit precision matrix
        X_sub_std = StandardScaler().fit_transform(X_sub)
        model = LedoitWolf()
        model.fit(X_sub_std)
        precision = model.precision_
        
        # Record edges above threshold
        for i in range(n_alleles):
            for j in range(i+1, n_alleles):
                if abs(precision[i, j]) > threshold:
                    edge_selections[(i, j)] += 1
    
    # Convert to stability scores (proportion of times selected)
    stability_scores = {edge: count / n_bootstrap for edge, count in edge_selections.items()}
    
    # Create stability score dataframe
    stability_data = []
    for (i, j), score in stability_scores.items():
        stability_data.append({
            'Allele_1': allele_list[i],
            'Allele_2': allele_list[j],
            'Family_1': family_list[i],
            'Family_2': family_list[j],
            'Stability_Score': score
        })
    
    stability_df = pd.DataFrame(stability_data)
    stability_df = stability_df.sort_values('Stability_Score', ascending=False)
    stability_df.to_csv("plots/edge_stability_scores.csv", index=False)
    
    print(f"\nStability selection complete.")
    print(f"Total edges detected at least once: {len(stability_scores)}")
    print(f"Edges with stability score >= 0.5: {sum(1 for s in stability_scores.values() if s >= 0.5)}")
    print(f"Edges with stability score >= 0.75: {sum(1 for s in stability_scores.values() if s >= 0.75)}")
    print(f"Edges with stability score >= 0.9: {sum(1 for s in stability_scores.values() if s >= 0.9)}")
    
    # Plot stability score distribution
    plt.figure(figsize=(8, 5))
    plt.hist(list(stability_scores.values()), bins=50, edgecolor='black')
    plt.axvline(x=0.5, color='r', linestyle='--', label='Score = 0.5', alpha=0.7)
    plt.axvline(x=0.75, color='orange', linestyle='--', label='Score = 0.75', alpha=0.7)
    plt.xlabel("Stability Score")
    plt.ylabel("Number of Edges")
    plt.title(f"Edge Stability Scores (n={n_bootstrap} bootstraps)")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig("plots/stability_score_distribution.pdf")
    plt.close()
    
    return stability_scores

def get_precision_matrix():
    X, allele_list, family_list = load_matrix()
    pX_std = StandardScaler().fit_transform(X)
    model = LedoitWolf()
    model.fit(pX_std)

    precision = model.precision_  # Inverse covariance matrix
    print ("Precision matrix shape:", precision.shape)
    print ("Precision matrix:", precision)

    # Generate sparsity path
    sparsity_df = estimate_cutoff(precision, allele_list)
    
    # Perform stability selection
    threshold = 0.02
    print(f"\n{'='*60}")
    print(f"STABILITY SELECTION ANALYSIS")
    print(f"{'='*60}")
    stability_scores = stability_selection(X, allele_list, family_list, 
                                          n_bootstrap=100, 
                                          subsample_ratio=0.8, 
                                          threshold=threshold)
    
    # Build network using both precision threshold and stability score
    print(f"\n{'='*60}")
    print(f"BUILDING NETWORK")
    print(f"{'='*60}")
    
    # Option 1: Use original threshold (0.02)
    threshold_edges = 0.02
    G = nx.Graph()
    node_family_list = {}
    
    for i in range(len(allele_list)):
        for j in range(i+1, len(allele_list)):
            if abs(precision[i, j]) > threshold_edges:
                # Check stability score if available
                edge_tuple = (i, j)
                stability = stability_scores.get(edge_tuple, 0.0)
                
                G.add_edge(allele_list[i], allele_list[j], 
                          weight=abs(precision[i, j]),
                          stability_score=stability)
                node_family_list[allele_list[i]] = family_list[i]
                node_family_list[allele_list[j]] = family_list[j]
    
    for node in node_family_list:
        G.add_node(node, label=node, type=node_family_list[node])
    
    nx.write_gml(G, "precision_matrix_network.gml")
    print(f"Network saved: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
    
    # Option 2: Build high-stability network (stability >= 0.75)
    G_stable = nx.Graph()
    node_family_list_stable = {}
    
    min_stability = 0.75
    for i in range(len(allele_list)):
        for j in range(i+1, len(allele_list)):
            edge_tuple = (i, j)
            stability = stability_scores.get(edge_tuple, 0.0)
            
            if stability >= min_stability:
                G_stable.add_edge(allele_list[i], allele_list[j], 
                                 weight=abs(precision[i, j]),
                                 stability_score=stability)
                node_family_list_stable[allele_list[i]] = family_list[i]
                node_family_list_stable[allele_list[j]] = family_list[j]
    
    for node in node_family_list_stable:
        G_stable.add_node(node, label=node, type=node_family_list_stable[node])
    
    nx.write_gml(G_stable, "precision_matrix_network_stable.gml")
    print(f"High-stability network saved (stability >= {min_stability}): "
          f"{G_stable.number_of_nodes()} nodes, {G_stable.number_of_edges()} edges")

# matrix = load_matrix()
# print (matrix)
get_precision_matrix()
# get_gene_precision_matrix()
