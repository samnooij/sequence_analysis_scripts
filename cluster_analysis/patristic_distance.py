#! /usr/bin/python

# Calculate patristic distances from tree file (in Newick format).

import dendropy

tree = dendropy.Tree.get(path = "tree.nwk", schema = "newick")

pdc = tree.phylogenetic_distance_matrix()

with open("tree.dist", 'w') as dist_matrix:
    dist_matrix.write("from;to;distance\n")
    for i, t1 in enumerate(tree.taxon_namespace[:-1]):
        for t2 in tree.taxon_namespace[i+1:]:
            dist_matrix.write("%s;%s;%s\n" % (t1.label, t2.label, pdc(t1, t2)))
