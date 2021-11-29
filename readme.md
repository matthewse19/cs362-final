# CS362: Computational Biology, Final Project
Matthew Smith-Erb
11/29/2021

This final project implements two phylogenetic consensus tree algorithms: Adams and Ml/majority/strict.

## Installation
Use the package manager [pip](https://pip.pypa.io/en/stable/) to install the following package [Biopython](https://biopython.org/)

    pip install biopython

## Files

 - consensus_cli.py - implements a CLI to run the consensus tree algorithms with input and output file
 - adams_tree.py - implements the Adams consensus tree algorithm
 - ml_tree.py - implements an Ml/majority/strict consensus tree algorithm
 - input_example.tree - an example file to use for input, with trees written on each line in [Newick](https://en.wikipedia.org/wiki/Newick_format) tree format

## consensus_cli.py Usage
usage: python consensus_cli.py [-h] [-d] [-l THRESHOLD] {Ml,adams} inputFile outputFile

Run Ml or Adams consensus tree algorithms

positional arguments:
  {Ml,adams} : Consensus algorithm to run to create the output tree
  inputFile : Name of file which includes the trees in Newick format
  outputFile : Name of the output file

optional arguments:
  -h, --help : show this help message and exit
  -d, --drawTree : Draw the output tree to the terminal with ascii characters
  -l THRESHOLD, --threshold THRESHOLD : When Ml algorithm is selected, this determines the ratio of input trees which must have the cluster, defaults to 0.5

### Examples
python consensus_cli.py  adams input_example.tree outputFile.tree

python consensus_cli.py -d Ml input_example.tree outputFile.tree

python consensus_cli.py -d -l=1 Ml input_example.tree outputFile.tree
