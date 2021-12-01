'''
Matthew Smith-Erb

11/29/2021

This file contains the implementation for the Ml consensus tree algorithm.
Note: bipartition is used here to refer to a specific cluster seperated from the rest of the taxa
'''

from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, Tree
from io import StringIO
from adams_tree import get_clade_names

def generate_ml_tree(trees, l):
    '''
    Takes in a list of trees, and a threshold l,
    returns the Ml tree.
    '''
    #only works for l in range [0.5, 1]
    assert (l >= 0.5 and l <= 1)

    #assert the taxa in each input tree is the same
    taxa_set = get_clade_names(trees[0].get_terminals())
    for tree in trees[1:]:
        assert (taxa_set == get_clade_names(tree.get_terminals())), "The input trees must have the same set of taxa!"

    #n number of taxa
    n = trees[0].count_terminals()

    #t number of trees
    t = len(trees)

    #dict where key is an immutable set of strings which correspond to a subset of the taxa
    #value is the ratio of input trees with the specific bipartition
    bipartition_count = {}
    #traverse each input tree post order, generating counts of each bipartition
    for tree in trees:
        for clade in tree.find_clades(order='postorder'):
            #create bipartition for the set of terminals that are descendants of this clade
            bipartition = set()
            if clade.is_terminal():
                bipartition = {clade.name}
            else:
                #create bipartition set from children's bipartition
                for child in clade:
                    child_bipartition = set(child.name.split(","))
                    bipartition = set.union(bipartition, child_bipartition)

                #store bipartition of internal node as the name of the internal node, to be used later
                clade.name = bipartition_to_string(bipartition)

            #make immutable so it can be used in dict
            bipartition = frozenset(bipartition)

            #increment occurances of this bipartition
            if bipartition in bipartition_count:
                bipartition_count[bipartition] += 1 / t
            else:
                bipartition_count[bipartition] = 1 / t

    #the keys represent bipartitions in the current Ml_tree, values are the parents of these nodes
    existing_bipartitions = {}
    #traverse each input tree in pre order, reconstructing the Ml from it
    for tree in trees:
        #c is the bipartition of the last node in input tree that corresponds to a majority/l bipartition that is an ancestor of the current node
        #initialize to be the bipartition at the root, so it includes all leaves
        c = frozenset(tree.root.name.split(","))
        for child in tree.root:
            preorder_reconstruct(child, c, bipartition_count, existing_bipartitions, l)

    #reverse the directions of the edges to make them go from root to descendents
    #key is bipartition/list of nodes and value is list of bipartitions that this key node points to
    Ml_dict = flip_paths(existing_bipartitions)

    #get the root of the tree
    root = find_root(existing_bipartitions)

    #converts the root and all descendents of the root to Clade objects
    converted_root = dict_to_clades(root, Ml_dict)
    Ml_tree = Tree(root = converted_root)
    relabel_tree(Ml_tree.find_clades())

    return Ml_tree

def preorder_reconstruct(clade, c, bipartition_count, existing_bipartitions, l):
    '''
    In preorder, check if the bipartition from the current clade is a majority bipartition,
    if it is a majority, add it to the existing_bipartitions or modify its parent if already an existing_bipartition.
    '''
    #get current bipartition from clade's name
    bipartition = frozenset(clade.name.split(","))
    count = bipartition_count[bipartition]

    #check if the bipartition is a majority
    if count > l or (count == l == 1):
        if bipartition in existing_bipartitions:
            #if the count of c is less than that of the parent in the Ml, switch the current clade's parent to be c

            parent = existing_bipartitions[bipartition]
            #TODOOOO seems fishy vvvv
            parent_count = bipartition_count[bipartition]
            #parent_count = bipartition_count[parent]

            c_count = bipartition_count[c]

            if parent_count > c_count:
                existing_bipartitions[bipartition] = c
        else:
            #bipartition doesn't exist in Ml tree, add it with its parent as c
            existing_bipartitions[bipartition] = c
        #set c to be the current clade's bipartition (which is a majority) and recursively pass it down
        c = bipartition

    #make recursive call to children
    for child in clade:
        preorder_reconstruct(child, c, bipartition_count, existing_bipartitions, l)

def flip_paths(graph):
    '''
    Given a set of edges from child to parent, return dict of edges from parent to children.
    Outputs a dict where key is the parent and value is list of children.
    '''
    corrected_tree = {}
    for out_neighbor in graph:
        in_neighbor = graph[out_neighbor]

        if in_neighbor in corrected_tree:
            corrected_tree[in_neighbor].append(out_neighbor)
        else:
            corrected_tree[in_neighbor] = [out_neighbor]

    return corrected_tree

def find_root(graph):
    '''
    Returns the root of the graph represented as a dict of edges from child (key) to parent (value)
    '''
    for in_neighbor in graph.values():
        if in_neighbor not in graph:
            return in_neighbor

def dict_to_clades(bipartition, dict):
    '''
    Converts a bipartition to a clade object and converts all of its descendents to clade objects
    '''
    clade_name = bipartition_to_string(bipartition)
    if bipartition in dict:
        #has children so convert them too
        return Clade(name=clade_name, clades = [dict_to_clades(child, dict) for child in dict[bipartition]])
    else:
        return Clade(name=clade_name)

def relabel_tree(clades):
    '''
    Set the length of the edges to zero and internal clades to have empty string names
    '''
    for clade in clades:
        if  not clade.is_terminal():
            clade.name = ""
        clade.branch_length = None

def bipartition_to_string(bipartition):
    '''
    Given a list of node names, concenatenate them with commas in between and return the string
    '''
    str = ""
    for terminal in bipartition:
        str += terminal + ","

    #remove last comma
    str = str[:-1]
    return str
