# St. Andrews: Machine Learning Phylogenomics


Welcome to this hands-on exercise in machine learning phylogenomics.<br>
In this assignment, you will train your own neural network to predict the topology of 4-taxon trees and 
then calculate the congruence between predicted quartets and a reference tree topology.



## Dataset and Scripts 

We will use a dataset [from this paper](https://doi.org/10.1186/1471-2148-11-114).<br>
The repository includes the following:<br>
Train_neural_network.py — script to train a neural network on site pattern frequencies.<br>
training_data/ — contains precomputed training data:<br>
>site_pattern_frequencies.npy: input features<br>
>topology_labels.npy: corresponding tree topologies<br>

concatenated_alignment_of_frogs_genes.fas — alignment of real frog gene sequences.<br>
frog_concat.part.contree — reference tree for quartet mapping.<br>
quartet_mapping.py — script to predict quartet topologies and compare them to the reference tree.<br>
visualise_tree.py — script to visualise a reference topology.<br>

## Training Dataset 

The training dataset was generated using  [PolyMoSim](https://github.com/cmayer/PolyMoSim/tree/main), a tool for simulating sequence evolution.<br> 
It includes 30,000 four-taxon alignments, each evolved under the GTR+I+G substitution model.<br>
The simulation parameters were chosen to broadly reflect the range of values observed in empirical datasets.

## Getting Started 
Let's start by cloning this repository:
  
```
git clone https://github.com/KulyaNikita/Machine-Learning-Phylogenomics-Workshop-St-Andrews.git
cd Machine-Learning-Phylogenomics-Workshop-St-Andrews
```

## Set Up the Environment

To install required packages, use:

```
mamba env create -f environment.yml

conda activate deepnnphylo 
```

To convert alignments into site pattern frequencies, compile the support script:
 
```
cd quartet_pattern_counter/
chmod u+x compile.sh
./compile.sh
cd ..
```
## Tree Visualization

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


## Predict and Map Quartets 

Use your trained model to predict quartet topologies and calculate quartets congruence with the reference tree:
```
python3 quartet_mapping.py frog_concat.part.contree concatenated_alignment_of_frogs_genes.fas MySuperCoolModel
```
Results will be saved in the results/ folder. 

What output files are generated?
What percentage of predicted quartets match the reference?
What does this say about the difficulty of the inferred tree?

Congratulations! Well done! 


