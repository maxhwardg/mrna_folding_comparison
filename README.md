# Comparison of mRNA Folding Packages
Code for the paper entitled "mRNA Folding Algorithms for Structure and Codon Optimization".

Includes Python scripts for running CDSfold, LinearDesign, and DERNA. This
comprises checking them for correctness, speed and memory usage. Although it is
not mentioned in the paper, this repository also has benchmarking code for
the mRNAfold package. Results are included in the data/ directory.
Code is located in src/, including bash scripts that can be used to generate
the plots in the paper.

Python dependencies can be found in requirements.txt. The code assumes the user
has downloaded and built the mRNA folding packages into a directory:
path_to_this_repo/extern. See the relevant package repository for build instructions.
Also, see src/bridge.py for specifics on how this code calls the executables for each
package.
