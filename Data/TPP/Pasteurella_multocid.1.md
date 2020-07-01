## TPP - Pasteuralla Multocid

Using the clustering pipeline, the CH index criteria leads to an edge 
sparsification threshold of 96%. A detailed comparison between the resulting 
clustering and the profiling suggests that while the two methods perform 
relatively comparably to capture the featured helix classes, community 
detection doesn't put the selected profiles separately into distinct clusters. 
Clustering is an improvement due to its retention of information from helix 
class 16 (which was not part of any selected profile). Helix class 16, 
which is incompatible with the featured helix class 87, appears in 7 of the 
34 clusters, which covers 12.4% of the sample. This example demonstrates 
that the clustering method can capture structural alternatives lost in 
profiling due to its threshold calculations for features and selected profiles.

This clustering was generated using an ensemble from RNAstructure with seed 1. 
It used the frequency distance metric between extended profiles and removed the 
lowest 96% of edges. GMC was used for clustering.