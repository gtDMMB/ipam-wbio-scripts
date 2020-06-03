import numpy as np
import score_helices

# Computes weighted symmetric difference of two structures
# first, second are structures (id, [helix classes])
# helices is a list of helix classes (id, (i,j,k), freq)
# weight should be one of "none", "freq", "length", "energy"
def Struct_Symdiff(first, second, helices, weight, sequence=None):

    # convert to sets and find classes in symm diff
    diff_set = set(first).symmetric_difference(set(second))

    # weighting "none": return length of the diff_list
    if weight == "none":
        return len(diff_set)
    # weight "freq": return symm_diff using weighted frequency
    elif weight == "freq":
        return sum(helices[idx]['freq'] for idx in diff_set)
    # weight "length": return symm_diff using weighted length
    elif weight == "length":
        return sum(helices[idx]['ijk'][2] for idx in diff_set)
    # weight "energy": return symm_diff using weighted energy
    elif weight == "energy" and sequence is not None:
        return sum(score_helices.helix_class_energy(
            [helices[idx] for idx in diff_set],sequence))
    print("Undefined weighting or no sequence provided; supported: {none, freq, length, energy}.")
    exit()

# Computes weighted symmetric differences for all pairs of structures
# Returns an array of size len(structures)*len(structures)
# structures and helices are (standard) lists of these structures
# weight should be one of "none", "freq", "length", "energy"
# sequence is required if weight is "energy"
def All_Struct_Symdiff(structures, helices, weight, sequence=None):
    #Create an empty 2D array
    symm_diff = np.empty((len(structures),len(structures)))

    for i in range(len(structures)):
        for j in range(len(structures)):
            if i == j:
                symm_diff[i][j] = 0
            elif i > j:
                symm_diff[j][i] = symm_diff[i][j] = Struct_Symdiff(
                    structures[i][1], structures[j][1], helices, weight, sequence)
    return symm_diff

if __name__ == "__main__":
    import sys
    import data
    import graph
    import networkx
    import helpers

    if len(sys.argv) != 4:
        print("Usage: python structure_similarity.py structure.out verbose.txt output_prefix")
        exit()

    structfile = sys.argv[1]
    verbosefile = sys.argv[2]
    file_prefix = sys.argv[3]

    structures = data.Read_Structures(structfile)
    helices = data.Read_HelixClasses(verbosefile)
    sequence = data.Read_Sequence(verbosefile)

    #Code for testing symmetric difference function

    symm_diff = All_Struct_Symdiff(structures, helices, "energy",sequence)
    similarity = helpers.Score_To_Sim(symm_diff)
    print("Symmetric Difference Matrix:")
    print(symm_diff)

    G = graph.mat_to_weighted_graph(
        similarity,max_edge_count=500,remove_isolated=True)
    graph.plot_graph_with_communities(
        G,
        show_node_labels=False,
        community_detection_func=networkx.algorithms.community.label_propagation.asyn_lpa_communities,
        filename=file_prefix + ".png")
    
