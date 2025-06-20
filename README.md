# St. Andrew: Machine Learning Phylogenomics


Welcome to this hands-on exercise in machine learning phylogenomics.
In this assignment, you will train your own neural network to predict the topology of 4-taxon trees and  
then calculate the congruence between predicted quartets and a reference tree topology.



## Objective and data


We will use a dataset [from this paper](https://doi.org/10.1186/1471-2148-11-114).
The repository includes the following:
Train_neural_network.py — script to train a neural network on site pattern frequencies.
training_data/ — contains precomputed training data:
    site_pattern_frequencies.npy: input features
    topology_labels.npy: corresponding tree topologies
concatenated_alignment_of_frogs_genes.fas — alignment of real frog gene sequences.
frog_concat.part.contree — reference tree for quartet mapping.
quartet_mapping.py — script to predict quartet topologies and compare them to the reference tree.
visualise_tree.py — script to visualise a reference topology

Let's start by cloning this repository:
  
```
git clone https://github.com/KulyaNikita/Machine-Learning-Phylogenomics-Workshop-St-Andrews.git
cd Machine-Learning-Phylogenomics-Workshop-St-Andrews
```

## Environment setting


Before starting a neural network training we should install all required packages. For that you should run commands below:

```
mamba env create -f environment.yml

conda activate deepnnphylo 
```

Also we need to activate a supportative script for a conversion of alignments to site pattern frequencies.
For that run the following commands
 
```
cd quartet_pattern_counter/
chmod u+x compile.sh
./compile.sh
cd ..
```
## Tree visualization

You can inspect the reference phylogeny using:

```
python3 visualise_tree.py frog_concat.part.contree

```
Are the branches long or short?

Now we are ready to start our analysis! 

## Neural Network Training 

Train the model using:
```
python3 Train_neural_network.py training_data/site_pattern_frequencies.npy training_data/topology_labels.npy > my_nn_training.txt 
```

What accuracy did your model achieve?


## Predict all the quartets of the alignment and map them back to the reference tree topology 

Use your trained model to predict quartet topologies and calculate quartets congruence with the reference tree:
```
python3 quartet_mapping.py frog_concat.part.contree concatenated_alignment_of_frogs_genes.fas MySuperCoolModel
```
Results will be saved in the results/ folder. 

What output files are generated?
What percentage of predicted quartets match the reference?
What does this say about the difficulty of the inferred tree?

Congratulations! Well done! 


