from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, Tree
from io import StringIO
import primes
import random

def generate_majority_tree(trees):
    #t number of trees
    t = len(trees)
    print("t/num of trees",t)

    #n number of taxa (assert for each tree)
    n = trees[0].count_terminals()
    for tree in trees[1:]:
        assert n == tree.count_terminals()
    taxa = [clade.name for clade in trees[0].get_terminals()]
    print("tax",taxa)
    print("n/num of taxa",n)

    #calculate m1
    m1 = primes.smallest_prime_greater_than(t * n)
    print("m1",m1)

    #calculate m2
    m2 = primes.primes[-1]
    print("m2", m2)

    #random list of n integers in range [0, m1) used for hash function
    a1 = random.sample(range(m1), k=n)
    a2 = random.sample(range(m2), k=n)

    hash_to_taxa = {}

    table = [[] for _ in range(m1)]
    tree_hashes = {}
    #traverse each input tree post order, generating counts of each bipartition
    for tree in trees:
        #TODO: check collisions
        postorder_nodes = tree.find_clades(order='postorder')
        hash_values = {}
        for clade in postorder_nodes:
            h1,h2 = 0,0
            if clade.is_terminal():
                leaf_id = taxa.index(clade.name)
                h1 = a1[leaf_id]
                h2 = a2[leaf_id]
                hash_to_taxa[(h1, h2)] = clade.name
                print(clade.name, h2)
                hash_values[clade] = (h1,h2)
            else:
                for child in clade:
                    hashes =  hash_values[child]
                    h1 += hashes[0]
                    h2 += hashes[1]
                h1 = h1 % m1
                h2 = h2 % m2
                hash_values[clade] = (h1,h2)
            bipartition_exists = False
            for entry in table[h1]:
                if entry[0] == h2:
                    entry[1] += 1
                    bipartition_exists = True
                    break
            if not(bipartition_exists):
                table[h1].append([h2, 1])
        tree_hashes[tree] = hash_values
    print(table)

    #traverse each input tree pre order, creating the Ml from it

    existing_bipartitions = {}
    for tree in trees:
        hash_values = tree_hashes[tree]
        c = hash_values[tree.root]
        for child in tree.root:
            preorder_reconstruct(child, c, hash_values, table, existing_bipartitions, t)
    print("Existing bipartitions", existing_bipartitions)
    Ml_dict = flip_paths(existing_bipartitions)
    root = find_root(existing_bipartitions)
    print("Ml_dict", Ml_dict)
    print(root)
    converted_root = dict_to_clades(root, Ml_dict)
    Ml_tree = Tree(root = converted_root)
    print(hash_to_taxa)
    relabel_tree(Ml_tree.find_clades(), hash_to_taxa)
    print(Ml_tree)
    Phylo.write(Ml_tree, 'idk.nhx', 'newick')
    Phylo.draw_ascii(Ml_tree)

def preorder_reconstruct(clade, c, hash_values, table, existing_bipartitions, t):
    hashes = hash_values[clade]
    count = count_from_hash(hashes, table)
    #TODO: generalize for l - RESTRICT TO BE BETWEEN 0.5 and 1 (count/t is greater than or EQUAL to this range)

    if count / t >= 1:
        if hashes in existing_bipartitions:
            #if the count of c is less than that of the parent, switch the current clade's parent to be c
            parent = existing_bipartitions[hashes]

            parent_count = count_from_hash(parent, table)

            c_count = count_from_hash(c, table)

            if parent_count > c_count:
                existing_bipartitions[hashes] = c
        else:
            #doesn't exist in Ml tree, add it with its parent as c
            existing_bipartitions[hashes] = c
        #set c to be the current clade's hashes and recursively pass it down
        c = hashes

    for child in clade:
        preorder_reconstruct(child, c, hash_values, table, existing_bipartitions, t)

def flip_paths(graph):
    #graph is a dict of paths in a tree going wrong way
    corrected_tree = {}
    for out_neighbor in graph:
        in_neighbor = graph[out_neighbor]

        if in_neighbor in corrected_tree:
            corrected_tree[in_neighbor].append(out_neighbor)
        else:
            corrected_tree[in_neighbor] = [out_neighbor]

    return corrected_tree

def find_root(graph):
    for in_neighbor in graph.values():
        if in_neighbor not in graph:
            return in_neighbor

def dict_to_clades(clade, dict):
    #TODO: make more self encompsing
    if clade in dict:
        #has children so convert them too
        return Clade(name=clade, clades = [dict_to_clades(child, dict) for child in dict[clade]])
    else:
        return Clade(name=clade)

def hash(input, m, a):
    sum = 0
    assert len(a) == len(input)
    for idx, a_i in enumerate(a):
        sum += (a_i * input[idx] % m)
    return sum

def count_from_hash(hashes, table):
    h1 = hashes[0]
    h2 = hashes[1]
    for entry in table[h1]:
        if entry[0] == h2:
            return entry[1]

def relabel_tree(clades, lookup):
    print(lookup)
    for clade in clades:
        if clade.is_terminal():
            clade.name = lookup[tuple(clade.name)]
        else:
            clade.name = ""
        clade.branch_length = None
if __name__ == "__main__":
    '''Phylo.read(StringIO("((A, (B, C)), (E, D))"), "newick"),
    Phylo.read(StringIO("((B, (A, C)), (E, D))"), "newick"),
    Phylo.read(StringIO("((C, (B, A)), (D, E))"), "newick")'''
    '''
    Phylo.read(StringIO("((1,(2,(3,4))),5)"), "newick"),
    Phylo.read(StringIO("(((1,2),(3,4)),5)"), "newick"),
    Phylo.read(StringIO("((((1,2),3),4),5)"), "newick")
    '''
    '''
    Phylo.read(StringIO("((A,B),(C,D))"), "newick"),
    Phylo.read(StringIO("((C,B),(A,D))"), "newick")
    '''
    trees = [
        Phylo.read(StringIO("((1,(2,(3,4))),5)"), "newick"),
        Phylo.read(StringIO("(((1,2),(3,4)),5)"), "newick"),
        Phylo.read(StringIO("((((1,2),3),4),5)"), "newick")

        #Phylo.read(StringIO("(S0, (S1, ((S2, S4), S3)))"), "newick")
    ]

    generate_majority_tree(trees)
