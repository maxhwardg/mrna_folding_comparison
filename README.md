# Comparison of mRNA Folding Packages
Code for the paper entitled "mRNA Folding Algorithms for Structure and Codon Optimization".

Relevant information from the paper:
"We conducted a series of experiments to compare existing mRNA folding software packages including LinearDesign, CDSfold, and DERNA. These were downloaded from their respective GitHub repositories using commits f0126ca, 06f3ee8, and ac84b6f compiled from source on Ubuntu 24.04 using GCC 13.2.0. All experiments were performed on the same Ubuntu system equipped with an AMD 7950X processor. The Homo sapiens codon frequency table from the Kazusa database was used for all experiments."


## Repository Structure
This repository contains Python scripts for running [CDSfold](https://github.com/gterai/CDSfold), [LinearDesign](https://github.com/LinearDesignSoftware/LinearDesign), and [DERNA](https://github.com/elkebir-group/derna).
This comprises checking them for correctness, speed and memory usage. Although it is not mentioned in the paper, this repository also has benchmarking code for the [mRNAfold](https://github.com/maxhwardg/mRNAfold) package.
Benchmark results are included in the data/ directory along with breaking cases the the codon frequency table.
Code is located in src/, including bash scripts that can be used to generate
the plots in the paper.

## Dependencies
Python dependencies can be found in requirements.txt.
The code assumes the user has downloaded and built the mRNA folding packages into a directory: path_to_this_repo/extern. This is not necessary for plotting existing results.
However, to run new benchmarks, the built packages need to be added.
See the relevant package repository for build instructions.
Also, see src/bridge.py for specifics on how this code calls the executables for each package.
