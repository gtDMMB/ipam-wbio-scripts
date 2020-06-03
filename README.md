
# Scripts for the 2019 IPAM Women in Biology Workshop Project in RNA Ensemble Analysis

These scripts are written for Python 3 and untested with 2.\*. They are known to 
work with Python 3.6, although they likely are functional with earlier versions. 
The Numpy, Scipy, Matplotlib, and NetworkX python packages are required.
They are dependent on two external software projects. Ensemble sampling is 
performed using the Mathews Lab 
[RNAstructure](https://rna.urmc.rochester.edu/RNAstructure.html) ensemble 
sampling component. The ensemble data is filtered and summarized using 
[RNAProfiling](http://rnaprofiling.gatech.edu/).

To run, use the ipam-fun.py script. For example
```python3 ipam-fun.py ~/tRNA/B0000234/B0000234_seq.txt ~/Desktop/RNAStructure/```

It is important that the input sequences are each in a folder with the name of the 
sequence, and the sequence files need to be of the form *_seq.txt. The script 
should work with either absolute or relative paths to the sequence file.
