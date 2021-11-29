'''
Matthew Smith-Erb

11/29/2021

This file contains the implementation for the Adams consensus tree algorithm
'''

import sys
from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, Tree
from io import StringIO
import copy

def generate_adams_tree(trees):
    '''
    Takes in a list of trees, asserts they have same taxa, and return Adams consensus tree.
    Assertion made in this helper so it isn't made every recursive call in the recursive_adams_tree() function
    '''
    taxa_set = get_clade_names(trees[0].get_terminals())
    for tree in trees[1:]:
        assert (taxa_set == get_clade_names(tree.get_terminals())), "The input trees must have the same set of taxa!"

    return recursive_adams_tree(trees)

def recursive_adams_tree(trees):
    '''
    Takes in a list of trees and returns the Adams consensus tree
    '''
    #base case: the trees only have 1 leaf
    if trees[0].count_terminals() == 1:
        return trees[0]


    partition_product = compute_partitions(trees)

    #list of trees which is the adams tree for the input trees under a restriction
    adams_trees = []

    for restriction in partition_product:
        #restrict each tree to the restriction
        restricted_trees = [restrict_tree(tree, restriction) for tree in trees]
        #get this adams tree for the restricted trees
        adams_tree = generate_adams_tree(restricted_trees)
        adams_trees.append(adams_tree)

    #create and return a tree whose root has a child for each of the restricted adams trees
    root = Clade(clades=[tree.root for tree in adams_trees])
    return Tree(root = root)

def compute_partitions(trees):
    '''
    Computes the partition product of the set of input trees.
    Returns a list of sets of partitions.
    '''
    leaves = trees[0].get_terminals()

    t = len(trees)

    partition_dict = {}

    #for each leaf, check which partion in each input tree it is in
    for leaf in leaves:
        #t length vector representing the leaf
        #leaf_vector[j] = i means the leaf is a descendent of the ith child of the root of input tree j.
        leaf_vector = [-1] * t
        for j, tree in enumerate(trees):
            for i, child in enumerate(tree.root):
                if leaf.name in get_clade_names(child.get_terminals()):
                    leaf_vector[j] = i
                    break

        #convert vector to tuple, so it can be used as a key
        leaf_tuple = tuple(leaf_vector)
        #add the leaf's name to a list at the leaf's tuple encoding
        if leaf_tuple in partition_dict:
            partition_dict[leaf_tuple].append(leaf.name)
        else:
            partition_dict[leaf_tuple] = [leaf.name]

    #return the set of leaf names at for each encoding
    return [set(b) for b in partition_dict.values()]

def restrict_tree(tree, restriction):
    '''
    Returns the given tree under the given restriction.
    This is the tree with only the paths to leaves in the restriction.
    This function paired with the get_parent() makes this Adams implementation no longer O(tn^2)
    '''

    tree = copy.deepcopy(tree)
    rescursive_restrict_tree(tree, restriction)

    #get rid of internal nodes which do not have branches
    root = tree.root
    for clade in tree.get_nonterminals():
        if len(clade) == 1 and clade is not root:
            #make the child of this clade be the child of this clade's parent
            child = clade.clades[0]
            parent = get_parent(tree, clade)
            parent.clades.remove(clade)
            parent.clades.append(child)

    #if tree has one leaf, return a tree with just a root, no edges
    if tree.count_terminals() == 1:
        return Tree(root = tree.get_terminals()[0])

    #collapse the one internal node that is the parent to all leaves
    if len(tree.root) == 1 and len(tree.root.clades[0]) != 1:
        tree.root.clades = tree.root.clades[0].clades

    return tree

def rescursive_restrict_tree(tree, restriction):
    '''
    Modifies the input tree to be restricted under the given restriction.
    Works recursively to removed unwanted nodes from the tree.
    '''
    #list of children of the root to later remove
    unwanted_children = []
    for child in tree.root:
        subtree_leaves = get_clade_names(Tree(root=child).get_terminals())

        #if the restriction set and the set of leaves of the subtree is disjoint, this child is unwanted
        if restriction.intersection(subtree_leaves) == set():
            unwanted_children.append(child)
        else:
            #otherwise, there is some overlap between subtree leaves and restriction, so recurse
            rescursive_restrict_tree(Tree(root=child), restriction)

    #remove these children from the root
    for unwanted in unwanted_children:
        tree.root.clades.remove(unwanted)

def get_parent(tree, child_clade):
    '''
    Given a tree and a child_clade, return the parent clade.
    Modified from:
        https://biopython.org/wiki/Phylo_cookbook
    Runs in lineare time with size of tree, making this Adams implementation not O(tn^2)
    '''
    node_path = tree.get_path(child_clade)
    #the path doesn't include the root by default, so if path is length 1, the parent is the root
    if len(node_path) == 1:
        return tree.root
    return node_path[-2]

def get_clade_names(clades):
    '''
    Returns a set of strings which are the names of the inputted list of clades
    '''
    return {clade.name for clade in clades}
