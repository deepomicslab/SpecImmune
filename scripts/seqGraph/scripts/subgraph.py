import networkx as nx
import sys
import os
def main(in_file, out_dir):
    fin = open(in_file)
    G = nx.Graph()
    for line in fin:
        values = line.split(" ")
        G.add_edge(values[0], values[2])
    sub_graphs = nx.connected_components(G)
    print(sub_graphs)

    i = 0
    for subg in sub_graphs:
        g = G.subgraph(subg)
        fout = open(os.path.join(out_dir,"graph"+str(i)), "w")
        for e in g.edges:
            fout.write(str(e))
        i = i+1
if __name__ == '__main__':
    main(sys.argv[1],sys.argv[2])