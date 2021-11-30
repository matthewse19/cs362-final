'''
Matthew Smith-Erb

11/29/2021

This scripts acts as a CLI to run the Adams and Majority consensus tree algorithms.
To get help on running the script, run 'python consensus_cli.py --help' or see the readme.md
'''

import argparse
from ml_tree import generate_ml_tree
from adams_tree import generate_adams_tree
from Bio.Phylo.BaseTree import Tree
from Bio import Phylo
from io import StringIO

def _restricted_float(x):
    '''
    Takes an input which would be a CLI input, returns it if it's a float and in [0.5,1]
    Modified from:
        https://stackoverflow.com/a/12117065
    '''
    try:
        x = float(x)
    except ValueError:
        raise argparse.ArgumentTypeError("%r not a floating-point literal" % (x,))

    if x < 0.5 or x > 1.0:
        raise argparse.ArgumentTypeError("%r not in range [0.5, 1.0]"%(x,))
    return x

def get_parsed_arguments():
    '''Creates and returns a parser object to get arguments:
    Algorithm name
    Input file name
    Output file name
    Optional flag for whether to print to terminal screen
    Optional threshold level for ml trees, defaults to 0.5 (a majority tree)
    '''

    parser = argparse.ArgumentParser(description='Run Ml or Adams consensus tree algorithms')
    parser.add_argument('algorithm', choices = ['Ml', 'adams'],
        help='Consensus algorithm to run to create the output tree')
    parser.add_argument('inputFile',
        help='Name of file which includes the trees in Newick format')
    parser.add_argument('outputFile',
        help='Name of the output file')
    parser.add_argument('-d', '--drawTree', action='store_true',
        help='Draw the output tree to the terminal with ascii characters')
    parser.add_argument('-l', '--threshold', type=_restricted_float, default = 0.5,
        help='When Ml algorithm is selected, this determines the minimum ratio of input trees which must have the cluster, defaults to 0.5')

    parsed_arguments = parser.parse_args()
    return parsed_arguments

if __name__ == '__main__':

    arguments = get_parsed_arguments()

    #get each line of the input file
    with open(arguments.inputFile) as file:
        input_lines = [line.rstrip() for line in file.readlines()]

    #read each line of the input file as a newick tree
    phylo_input_trees = []
    for line in input_lines:
        tree = Phylo.read(StringIO(line), 'newick')
        phylo_input_trees.append(tree)


    #run the correct algorithm
    if arguments.algorithm == 'adams':
        output_tree = generate_adams_tree(phylo_input_trees)
    elif arguments.algorithm == 'Ml':
        output_tree = generate_ml_tree(phylo_input_trees, arguments.threshold)

    Phylo.write(output_tree, arguments.outputFile, 'newick')

    #print tree if flag is present
    if arguments.drawTree:
        print("Consensus algorithm: " + arguments.algorithm)
        Phylo.draw_ascii(output_tree)
