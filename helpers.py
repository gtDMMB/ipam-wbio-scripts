import numpy as np
import matplotlib.pyplot as plt
from collections import Counter

#Sort a list by the specified tuple element (tup_index).
#Does not modify original list.
def Sort_Tuples(list, tup_index):
    return(sorted(list, key = lambda tup: tup[tup_index]))

#Find max structure size
def max_hc(structures):
    return len((max(structures, key=lambda x:len(x[1]))[1]))

#iterates over all of the base pairs in a helix class
def helix_class_iterator(helix_class, sequence,sort=False):
    if not isinstance(helix_class['ijk'],tuple):
        print("Bad input to helix_class_iterator probably")
    i,j,k = helix_class['ijk']
    for idx in range(k):
        #offset by one because indexing
        if sort:
            yield tuple(sorted([sequence[i + idx - 1], sequence[j - idx - 1]]))
        else:
            yield (sequence[i + idx - 1], sequence[j - idx - 1])

#a function for taking a list of structures and removing duplicates.
#it returns a list of tuples where the first element in each tuple is the count
#and the second element is the structure. This should allow it to be
#interoperable with most of the code which is already written.
#it also works for deduplicating profiles
def Deduplicate_structures(structures):
    counted = Counter([tuple(sorted(struct[1])) for struct in structures])
    out = [(value,key) for key, value in counted.items()]
    return sorted(out,reverse=True)

#Plot histogram of entries in symmetric 2D score_array (e.g. of pairwise similarities)
#Currently ignores main diagonal (assumes 0)
#Second parameter is number of bins
def Array_Histogram(score_array, bins, show=False, save=True, filename='plot.png'):
    #Create 1-d list to store similarities for histogram generation
    list_triu = []
    for i in range(len(score_array)-1):
	       list_triu.extend(score_array[i, i+1:len(score_array)])
    #Create plot
    n, bins, patches = plt.hist(list_triu, bins)
    if(show):
        plt.show()
    if(save):
        plt.savefig(filename)
    plt.close()

#calculate entry-wise reciprocal of array (set 0's to 1)
def Score_To_Sim(score_array):
    sim = np.copy(score_array)
    mask = sim !=0 #choose non-zeros
    sim[mask] = 1./sim[mask] #take reciprocals
    mask = sim == 0 #choose zeros
    sim[mask] = -1 #set to -1
    return sim

#function to partition a graph using a given networkx community detection function
def cluster_graph(G,community_detection_func,edge_weight_key='weight'):
    try:
        return community_detection_func(G,weight=edge_weight_key)
    except ZeroDivisionError as e:
        return [[node] for node in G] #There was an issue with the community detection

if __name__ == "__main__":
    import sys
    import data
    import graph
    import structure_similarity

    if len(sys.argv) != 4:
        print("Usage: python helpers.py structure.out verbose.txt output_prefix")
        exit()

    structfile = sys.argv[1]
    verbosefile = sys.argv[2]
    file_prefix = sys.argv[3]

    structures = data.Read_Structures(structfile)
    helices = data.Read_HelixClasses(verbosefile)

    structures_set = Deduplicate_structures(structures)

    #Test max structure length
    print("Max structure length:")
    print(max_hc(structures))

    print("Deduplicated structures:")
    print(structures_set)

    #Test array histogram
    num_bins = 50
    symm_diff = structure_similarity.All_Struct_Symdiff(structures, helices, "none")
    Array_Histogram(symm_diff, num_bins, False, True, file_prefix + '_hist_none.png')

    symm_diff = structure_similarity.All_Struct_Symdiff(structures, helices, "freq")
    Array_Histogram(symm_diff, num_bins, False, True, file_prefix + '_hist_freq.png')

    symm_diff = structure_similarity.All_Struct_Symdiff(structures, helices, "length")
    Array_Histogram(symm_diff, num_bins,False, True, file_prefix + '_hist_length.png')
