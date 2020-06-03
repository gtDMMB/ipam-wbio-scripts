import sys
import os

import networkx as nx
import matplotlib.pyplot as plt
from networkx.algorithms import bipartite
from networkx.algorithms import community
from networkx import edge_betweenness_centrality as betweenness
import numpy as np

#Import functionality
import data
import helpers
import structure_similarity
import score_helices
import graph
import cluster

#Check arguments to get structures and verbose output files
if len(sys.argv) < 2:
    print("IPAM structure similarity and clustering experiments")
    print("Usage: python ipam_fun.py path/to/sequence_file.txt [path/to/RNAStructure/directory/] [seed_file.txt]")
    print("An example of an RNAstructure directory is ~/Desktop/command-line/RNAstructure")
    print("The RNAStructure directory must be provided if samples have not already been generated")
    print("    in the proper directory structure")
    print("The RNAStructure seed is the index of the sample. So the seed for sample 10 is 10.")
    print("This script must be run from the directory it is in, otherwise it will not work")
    print("seed_file.txt is a series of rows of seeds. Each row should contain a single integer seed.")
    exit()

general_file_prefix = data.Get_File_Prefix(sys.argv[1])
output_directory = "./test_files/{}/".format(general_file_prefix)

rnastruct_dir = ""
seed_dir = ""
if len(sys.argv) > 2 and os.path.isdir(sys.argv[2]):
    rnastruct_dir = sys.argv[2]
    if len(sys.argv) > 3:
        seed_dir = sys.argv[3]
elif len(sys.argv) > 3 and os.path.isdir(sys.argv[3]):
    rnastruct_dir = sys.argv[3]
    seed_dir = sys.argv[2]
if rnastruct_dir == "":
    print("Warning: no RNAStructure directory was provided. This could cause a crash if new seeds or sequences are tested.")

num_bins=20
independent_samples = 1
difference_metrics = ["none","freq","length","energy"]
test_percentiles = [10, 20]

#Defining community detection functions to be used later

def community_gmc(G,weight_key='weight',**kwargs):
    return community.greedy_modularity_communities(G,weight=weight_key)

def community_lpa(G,**kwargs):
   return list(community.label_propagation_communities(G))

#will only find two communities (can be modified to find k communitites for a fixed k)
def community_gn(G,weight_key='weight',**kwargs):
    def most_central_edge(G):
        centrality = betweenness(G, weight=weight_key)
        return max(centrality, key=centrality.get)

    girvan_results = community.girvan_newman(G,most_valuable_edge=most_central_edge)
    return next(girvan_results)

community_functions = [(community_gmc, "gmc"),
                        (community_gn, "gn"),
                        (community_lpa,"lpa")]

if seed_dir == "":
    seed_list = [x+1 for x in range(independent_samples)]
else:
    with open(seed_dir,'r') as f:
        lines = f.read().splitlines()
    seed_list = [int(x) for x in lines]

#Iterating over samples, creates an independent set of structures for a single sequence each time
for iteration,seed in enumerate(seed_list):
    print("Starting sample number {} with seed {}".format(iteration,seed))

    file_prefix = general_file_prefix + "_{}".format(seed)

    #Uses the sample index as the seed for RNAstructure
    sequence, helices, structures, profiles =\
        data.Load_All_Data(sys.argv[1],rnastruct_dir,
		recalculate_sample=False,seed=seed)

    print("Number of unique structures: {}".format(len(structures)))

    #Store multiplicities in a list
    multiplicity = [x[0] for x in structures]

    #Iterate over the distance metrics for weighting symmetric difference
    for metric in difference_metrics:
        print("Starting calculations for metric {}".format(metric))

        #Calculate symmetric difference score
        symdiff = structure_similarity.All_Struct_Symdiff(structures, helices, metric, sequence)
        #Convert to similarity
        similarity = helpers.Score_To_Sim(symdiff)
        helpers.Array_Histogram(similarity, num_bins, False, True, 
            output_directory + file_prefix + '_hist_{}.png'.format(metric))

        #Iterate over percentiles and community detection functions
        for percentile in test_percentiles:
            for cluster_func, func_name in community_functions:
                print("Clustering {} percent of edges with {}".format(percentile,func_name))

                #Reasonable thresholds: fixed number (percentage) of edges
                #top X percent of similarities.
                G = graph.mat_to_weighted_graph(
                    similarity,threshold=np.percentile(similarity,percentile),
                    max_edge_count=5000,remove_isolated=False)

                partition = helpers.cluster_graph(G,cluster_func)
                ch_index = cluster.calculate_ch_index(structures,partition,helices,
                        sequence,verbose=False,metric=metric)

                print("{} cluster(s) found with a ch-index of {}".format(len(partition),ch_index))

                data.Save_Partition(
                    output_directory + file_prefix+"-{}-{}-{}-partition.txt".format(metric,percentile,func_name),
                    partition,data=structures,sequence=sequence,ch_index=ch_index)

                graph.plot_graph_with_communities(
                    G,
                    filename=output_directory + file_prefix+\
                        '-{}-{}-{}-graph.png'.format(metric,percentile,func_name),
                    node_labels=multiplicity,
                    node_partition=partition)

