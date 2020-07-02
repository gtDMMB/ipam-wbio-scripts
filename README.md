
# Scripts and Data for the 2019 IPAM Women in Biology Workshop Project in RNA Ensemble Analysis

## Description

This repository contains the scripts created for the 2019 IPAM Women in Biology 
Workshop Project in RNA Ensemble Analysis as well as the resulting data. 
`ipam-fun.py` is the main script, 
although many of the other files are runnable. Graph based clustering is performed 
on filtered sampled RNA secondary structural ensembles to identify important 
clusters of features. The output for various sequences as well as the sampled 
ensembles are in the `Data` folder. For more details, please see the associated 
book chapter (link will be made available).

## Dependencies

These scripts are written for Python 3 and untested with 2.\*. They are known to 
work with Python 3.6, although they likely are functional with earlier versions. 
The [Numpy](https://numpy.org/), [Scipy](https://www.scipy.org/), 
[Matplotlib](https://matplotlib.org/), and [NetworkX](https://networkx.github.io/) 
python packages are required. The scripts are dependent on two other software 
projects. Ensemble sampling is performed using the Mathews Lab 
[RNAstructure](https://rna.urmc.rochester.edu/RNAstructure.html) ensemble 
sampling component. The ensemble data is filtered and summarized using 
[RNAProfiling](http://rnaprofiling.gatech.edu/). The compiled RNAProfiling executable 
must be placed in a parallel directory to ipam-wbio-scripts. If this git repository 
was cloned to `~/Documents/replication/ipam-wbio-scripts` then there must be an 
executable `~/Documents/replication/src/RNAprofile`.

## Usage

To run, use the ipam-fun.py script. For example, if B0000234_seq.txt is a file containing 
a tRNA sequence and RNAStructure is placed on the user's desktop:

```python3 ipam-fun.py ~/tRNA/B0000234/B0000234_seq.txt ~/Desktop/RNAStructure/```

It is important that the sequence files be of the form *\_seq.txt. The script 
should work with either absolute or relative paths to the sequence file. An optional 
third argument is a file which contains line separated integer seeds. If this is not 
provided, then the number 1 is used as a single seed.

## Data 

For each sequence used in the final analysis, we provide the sequence, 
sampled ensemble, RNAProfiling output, visualization of the clusters, a text file containing the 
clusters, and a file containing a  description of the clusters and experimental parameters used to 
run the scripts.