## TPP - Thermoplasma Acidophilum

In the clustering pipeline, the CH index criteria monotonically increases as 
we sparsify, so the edge sparsification threshold is set to 99%. The 
clustering does not recapitulate the profiling results since the selected 
profiles are not separated into distinct clusters. The resulting clustering 
partition not only captures the same strong signal as profiling, but also 
retains information from features 2,12,17,32,74. By carefully analyzing the 
clustering, we also observed that helix class 43 is incompatible with helix 
classes 12,74, while helix class 18 is incompatible with helix classes 
17,32. Similarly to Pasteurella Multocid, the clustering of this example also 
retains structural alternatives lost in profiling.

This clustering was generated using an ensemble from RNAstructure with seed 1. 
It used the frequency distance metric between extended profiles and removed the 
lowest 99% of edges. GMC was used for clustering.