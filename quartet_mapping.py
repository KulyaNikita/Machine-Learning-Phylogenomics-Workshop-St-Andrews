#!/usr/bin/env python3
import sys
from ete4 import PhyloTree
from ete4.smartview import Layout
import re
import os
import subprocess
import numpy as np
from keras.models import load_model
from numpy import argmax
import itertools
from math import comb
from Bio import SeqIO


#Paths to the files and other arguments
ref_topology_path = sys.argv[1]
alignment_path = sys.argv[2]
machine_learning_model_path = sys.argv[3]

# Some helpful functions
# Some functions for species names extractions
def get_headers_from_indices(taxa_list,indices_tuple):
    # Get names from the list based on indices by subtracting 1
    taxa_from_indices = []
    for i in indices_tuple:
        taxa_from_indices.append(taxa_list[int(i)-1])

    return taxa_from_indices

def tree_topology(topology, header_list):
    if topology == 0:
        #        header_list = read_alignment_file(alignment)
        string = (
            "(("
            + str(header_list[0])
            + ","
            + str(header_list[1])
            + "),("
            + str(header_list[2])
            + ","
            + str(header_list[3])
            + "));"
        )

        return string
    elif topology == 1:
        #        header_list = read_alignment_file(alignment)
        string = (
            "(("
            + str(header_list[0])
            + ","
            + str(header_list[2])
            + "),("
            + str(header_list[1])
            + ","
            + str(header_list[3])
            + "));"
        )

        return string
    elif topology == 2:
        #        header_list = read_alignment_file(alignment)
        string = (
            "(("
            + str(header_list[0])
            + ","
            + str(header_list[3])
            + "),("
            + str(header_list[2])
            + ","
            + str(header_list[1])
            + "));"
        )

        return string

def quartet_mapping_nn(f,t_ref):
    failure_counts = 0
    success_counts = 0

    for quartet in f:
        t_quar = PhyloTree(quartet)
        list_of_quartet_sps = []
        for leaf in t_quar.leaves():
            list_of_quartet_sps.append(leaf.name)

        t_modified = t_ref.copy()

        t_modified.prune(list_of_quartet_sps)
        rf_dist = t_quar.compare(t_modified, unrooted=True)

        rf_value = rf_dist["rf"]

        if int(rf_value) == 0:
            success_counts += 1
        else:
            failure_counts += 1

    return success_counts, failure_counts

# MAKE OUTPUT DIRECTORY
output_dir = "results"
os.makedirs(output_dir, exist_ok=True)

# Create numpy array with site pattern frequencies
numpy_file = os.path.splitext(os.path.basename(alignment_path))[0]
numpy_output_file = os.path.join(output_dir, f"{numpy_file}_site_pattern_frequencies.npy")

command = f"quartet_pattern_counter/quartet-pattern-counter-v1.4  {alignment_path} {numpy_output_file}"
subprocess.run(command, shell=True)

# Load the memory-mapped array
frequency_array = np.load(numpy_output_file, mmap_mode='r')

# Load Machine Learning Model
machine_learning_model = load_model(machine_learning_model_path)

# Make Predictions
topology_prediction = machine_learning_model.predict(frequency_array)
list_predictions = argmax(topology_prediction, axis=1)

# Find number of sequences in the alignment
headers = [record.id for record in SeqIO.parse(alignment_path, "fasta")]

# all indices
all_indices = []
for combination in itertools.combinations(range(1, len(headers) + 1), 4):
    all_indices.append(combination)
# ALL QUARTETS FOR MAPPING
f = open("results/AllMappingQuartets.tre", "w")
# ITERATE THROUGH HEADER's indices AND predicted topologies and create a file with all tree topologies
for top_pred, alignment_indices in zip(list_predictions, all_indices):
    quartet_headers = get_headers_from_indices(headers, alignment_indices)
    one_topo = tree_topology(top_pred, quartet_headers)
    f.write(one_topo)
    f.write("\n")
f.close()


# Load ref phylogenetic tree
t_ref = PhyloTree(ref_topology_path)
f = open("results/AllMappingQuartets.tre", "r")
number_of_corr, number_of_fails = quartet_mapping_nn(f, t_ref)
f.close()
total_combinations = comb(len(headers), 4)
percentage_correct = (number_of_corr / total_combinations) * 100
percentage_incorrect = 100 - percentage_correct

with open("results/nn_mapping_results.txt", "a") as result_file:
    result_file.write(f"Percentage of the matched quartets is: {percentage_correct:.2f}%\n")
    result_file.write(f"Percentage of not matched quartets is: {percentage_incorrect:.2f}%\n")

