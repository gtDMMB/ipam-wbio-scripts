import networkx as nx
import matplotlib.pyplot as plt
from networkx.algorithms import bipartite
import numpy as np

#Returns a bipartite graph of helix classes versus structures, where edges represent containment
# Only considers the first num_structs elements of structures list.
def hc_structure_bipgraph(structures, num_structs):
    B = nx.Graph()
    max_structures=min(num_structs-1, len(structures))

    relevant_helices = set(
        helix for struct in structures[:max_structures] for helix in struct[1])

    edges = []
    #add nodes for each helix class and each structure
    for helix in relevant_helices:
        B.add_node("h"+str(helix),bipartite=0)
    for item in structures[0:max_structures]: #currently only first max_structures
        #add edges from each structure node to the helix nodes the structure contains
        B.add_node("s"+str(item[0]),bipartite=1)
        edges += [("h" + str(helix),"s"+str(item[0])) for helix in item[1]]
    B.add_edges_from(edges)
    return B

#Plots graph assuming that it is bipartite
def plot_bipgraph(graph):
    top_nodes = set(n for n,d in graph.nodes(data='bipartite') if d==0)
    color_list = ['#aa4444','#8866aa']
    colors = [color_list[node[1]] for node in graph.nodes(data='bipartite')]
    pos = nx.bipartite_layout(graph,top_nodes)
    nx.draw_networkx(graph,pos,with_labels=True,node_size=500,node_color=colors)
    plt.axis('off')
    plt.show()

#Create edge-weighted graph from symmetric matrix
# Setting threshold keeps only edges with similarity at least this value.
#if plot is set to true, we visualize with matplotlib
def mat_to_weighted_graph(
        score_array,
        threshold=0, #threshold for minimum edge weights
        max_edge_count=5000, #max number of edges to keep, may keep fewer to avoid arbitrary decisions
        remove_isolated=False): #removes vertices with no edges

    if np.count_nonzero(score_array) > 2*max_edge_count:
        threshold = max(
            threshold, np.partition(score_array.flatten(), -max_edge_count)[-max_edge_count])
    score_array = (score_array > threshold) * score_array

    G=nx.from_numpy_array(score_array)

    if remove_isolated:
        G.remove_nodes_from(list(nx.isolates(G)))
        G = nx.convert_node_labels_to_integers(G)

    return G

#plot a networkx graph and color nodes based on communities or input partitions
def plot_graph_with_communities(
        G, #graph to plot
        filename=None,#where to save the plot image. If blank, does not save.
        edge_weight_key='weight', #key for edge attribute to use as label (probably 'weight')
        node_labels=None, #labels for nodes as list (none if left empty)
        show_node_labels=True, #ignored if node_labels are provided
        plot=False, #whether to visualize the result
        community_detection_func=None, #a networkx community detection function
        node_partition=None): #a list of lists of node indices for partitions.
                              #node_partion is ignored when community_detection_func is provided

    node_colors = np.zeros(G.order()).astype(int)
    if community_detection_func is not None:
        communities = community_detection_func(G,weight=edge_weight_key)
        for i, community in enumerate(communities):
            node_colors[np.array(list(community)).astype(int)] = i + 1;
    elif node_partition is not None:
        for i, partition in enumerate(node_partition):
            node_colors[np.array(list(partition)).astype(int)] = i + 1
    if node_labels is not None:
        show_node_labels = False

    pos=nx.spring_layout(G,k=.5)
    nx.draw_networkx(G,pos,with_labels=show_node_labels,node_color=node_colors,cmap=plt.cm.Set1)

    if(node_labels is not None):
        labels = dict(enumerate(node_labels))
        nx.draw_networkx_labels(G, pos,labels=labels,font_size=8)

    plt.axis('off')
    if filename is not None:
        plt.savefig(filename,dpi=300)
    if filename is None or plot == True:
        plt.show()
    plt.close()

if __name__ == "__main__":
    import sys
    import data
    import structure_similarity
    import helpers

    #Check arguments to get structures, verbose output file, prefix for writing
    if len(sys.argv) != 4:
        print("Usage: python graph.py structures.out verbose.txt output_prefix")
        exit()

    #Read data from files
    structFile = sys.argv[1]
    verboseFile = sys.argv[2]
    file_prefix = sys.argv[3]
    structures = data.Read_Structures(structFile)
    helices = data.Read_HelixClasses(verboseFile)

    #B = hc_structure_bipgraph(structures,10)
    #plot_bipgraph(B)

    #Warning, this is slow
    symm_diff = structure_similarity.All_Struct_Symdiff(structures, helices, "none")
    similarity = (-1)**symm_diff + helpers.max_hc(structures)
    G = mat_to_weighted_graph(similarity, True) # not interesting right now
    plot_graph_with_communities(G)
