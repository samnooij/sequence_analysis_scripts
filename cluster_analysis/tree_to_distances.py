#! /usr/bin/python

#A script to calculate a pairwise distance matrix from a newick tree file.
# Made by Will: https://www.biostars.org/p/55410/#55416
# and slightly modified to correct typos.

from Bio import Phylo
from itertools import combinations

tree = Phylo.read("tree-file.nwk"), 'newick')
seq_names = tree.get_terminals()
dmat = {}
for p1, p2 in combinations(seq_names, 2):
    d = tree.distance(p1, p2)
    dmat[p1.name, p2.name] = d
    dmat[p2.name, p1.name] = d
