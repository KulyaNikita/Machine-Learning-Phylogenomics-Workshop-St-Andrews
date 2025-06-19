#!/usr/bin/env python3
from ete4 import Tree
import sys

'''
usage: python3 visualise_newick.py  newick_tree
'''
tree_file = sys.argv[1]

tree = Tree(tree_file, parser=1)

tree.explore()

input()
