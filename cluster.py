from structure_similarity import All_Struct_Symdiff, Struct_Symdiff
from random import randint
import sys

# COMPUTE CH-INDEX FOR A GIVEN METRIC
# Note: ensemble/cluster centroid is internally stored as (multiplicity, list helix classes)

def findCentroid(structures,helices,metric,sequence,verbose):
    sym_diff = All_Struct_Symdiff(structures,helices,metric,sequence)

    struct_scores = []
    for struct, differences in zip(structures, sym_diff):
        #get total distance to all other structures with weighting from counts
        total = sum(item[0]*diff for item, diff in zip(structures,differences))
        struct_scores.append((total,struct[1]))
    min_score = min(struct_scores)

    #find the number of unique structures with the minimum score
    centroid_count = sum(min_score[0]==score[0] for score in struct_scores)
    if (centroid_count > 1 and verbose):
        print("Different structures realize the minimum distance in the ensemble")

    return min_score[1]

def get_sorted_clusters(structures,partition):
    clusters = [[structures[idx] for idx in struct_list] for struct_list in partition]
    cluster_lengths = [sum(struct[0] for struct in cluster) for cluster in clusters]

    cluster_lengths, clusters = zip(*[(length, cluster) for length, cluster in\
            reversed(sorted(zip(cluster_lengths, clusters)))])

    return clusters, cluster_lengths

def calculate_ch_index(structures,partition,helices,sequence,verbose=False,metric="none"):
    ens_centroid = findCentroid(structures,helices,metric,sequence,verbose)

    clusters, cluster_lengths = get_sorted_clusters(structures,partition)
    cluster_centroids = [findCentroid(cluster,helices,metric,sequence,verbose) for cluster in clusters]

    B = 0
    W = 0
    for cluster, centroid, length in zip(clusters,cluster_centroids,cluster_lengths):
        add = length * Struct_Symdiff(ens_centroid,centroid,helices,metric,sequence)
        B += add

        W += sum(struct[0]*Struct_Symdiff(centroid,struct[1],helices,metric,sequence) 
                for struct in cluster)

    #struct[0] contains the number of repetitions of the structure
    n = sum(struct[0] for struct in structures)
    k = len(clusters)

    if (n == k):
        if verbose: print("Only one structure per cluster")
        return 0
    if (W == 0):
        if verbose: print("Zero within cluster variation")
        return 0
    if (k == 1):
        if verbose: print("Only one cluster provided")
        return 0

    ch_index = B / (k - 1) / (W / (n - k))
    return ch_index

if __name__ == "__main__":
    if (len(sys.argv) < 3):
        print("Usage: python ch_index path/to/sequence_file path/to/rnastructure [seed]")
        exit()

    sequence_file = sys.argv[1]
    rnastruct_dir = sys.argv[2]

    seed = 1
    if (len(sys.argv) == 4):
        seed = int(sys.argv[3])

    from data import Load_All_Data
    #Uses the sample index as the seed for RNAstructure
    sequence, helices, structures, profiles =\
        Load_All_Data(sequence_file,rnastruct_dir,
                    recalculate_sample=False,seed=seed)

    cluster_count = 5
    partition = [[] for _ in range(cluster_count)]
    [partition[randint(0,cluster_count-1)].append(i) for i in range(len(structures))]

    ch_index = calculate_ch_index(structures,partition,helices,sequence)

    print(ch_index)
    
    '''
    sequence = ['u']*10
    helices = [{'id':0,'ijk':(0,5,2),'freq':10},
            {'id':1,'ijk':(2,8,3),'freq':5},
            {'id':2,'ijk':(4,6,1),'freq':1},
            {'id':3,'ijk':(5,9,2),'freq':3}]

    structures = [(1,[0,1,2]),(1,[1,2,3]),(3,[2,3])]
    partition = [[0],[1,2]]

    ch_index = calculate_ch_index(structures,partition,helices,sequence)
    print(ch_index)
    '''

