# St. Andrew: Machine Learning Phylogenomics


This little exercise will introduce you to the basics of machine learning phylogenomics.

You will train your owm neural network model, allowing you to predict topologies of 4-taxa trees. 

Using this model you will estimate the difficultness of the provided tree topology.



## Objective and data


We will use a dataset [from this paper](https://doi.org/10.1186/1471-2148-11-114). In the repo we have following files (1) a script to train a neural network "Train_neural_network.py"
(2) beforehand prepared training data with saved site pattern frequencies and corresponding tree topologies under the folder "training_data/".
(3) the alignment file "concatenated_alignment_of_frogs_genes.fas"
(4) reference tree topology against which we will run quartet mapping "frog_concat.part.contree"
(5) a script to make predictions of the quartets and map them on the reference topology "quartet_mapping.py"

Let's start by cloning this repository:
  
```
git clone https://github.com/KulyaNikita/Machine-Learning-Phylogenomics-Workshop-St-Andrews.git

```

## Environment setting


Before starting a neural network training we should install all required packages. For that you should run commands below:

```
mamba env create -f environment.yml

conda activate deepnnphylo 
```

Also we need to activate a supportative script for a conversion of alignments to site pattern frequencies.
For that run the following commands

<details>
  <summary>Need help?</summary>
  
```
cd quartet_pattern_counter/
chmod u+x compile.sh
./compile.sh
```
</details>


Let's fix sequence names to get tidy files and trees! Also, having homogeneous names across ortholog groups is ESSENTIAL for the concatenation step. You can see that headers have the following format: `Genus_species_GENE_XXXX`. Can you simplify these with some `bash` commands so that only the common part is retained (i.e. `Genus_species`)?

<details>
  <summary>Need help?</summary>
  
```
for f in *fa; do awk -F"_" '/>/ {print $1"_"$(NF-2)}; !/>/ {print $0}' $f> out; mv out $f; done

# Explanation: loop through the files, modify the input file and save it to "out", overwrite the input with it. 
# Awk: using the "_" field separator, modify lines starting with ">" (sequence names) so that they will only 
# contain the second and second last elements separated by an underscore (>Genus_species).
# Print all other lines not starting with ">" (i.e. sequences) without modification.

# An easier version:
for f in *fa; do sed -E '/>/ s/_GENE.+//g' $f > out; mv out $f; done
# Explanation: in headers (lines starting with ">"), remove everything after "_GENE"
# Can you guess the difference with using awk? Look at Canis_lupus_familiaris
```
</details>

Don't forget to check the output: is your command doing what you want?


**NOTE ABOUT ORTHOLOGY**: Ensuring orthology is a difficult issue and often using a tool like Orthofinder might not be enough. Paralogy is tricky business! Research has shown (e.g. [here](https://www.nature.com/articles/s41559-017-0126) or [here](https://academic.oup.com/mbe/article/36/6/1344/5418531)) that including paralogs can bias the phylogenetic relationship and molecular clock estimates, particularly when phylogenetic signal is weak. Paralogs should always be removed prior to phylogenetic inference. But identifying them can be difficult and time consuming. One could build single-gene trees and look for sequences producing extremely long branches or clustering outside of the remaining sequences.



## Pre-alignment quality filtering


Often, transcriptomes and genomes have stretches of erroneous, non-homologous amino acids or nucleotides, produced by sequencing errors, assembly errors, or errors in genome annotation. But until recently, these type of errors had been mostly ignored because [no automatic tool could deal with them](https://natureecoevocommunity.nature.com/users/54859-iker-irisarri/posts/37479-automated-removal-of-non-homologous-sequence-stretches-in-phylogenomic-datasets).

We will use [PREQUAL](https://doi.org/10.1093/bioinformatics/bty448), a new software that takes sets of (homologous) unaligned sequences and identifies sequence stretches sharing no evidence of (residue) homology, which are then masked in the output. Note that homology can be invoked at the level of sequences as well as of residues (amino acids or nucleotides). Running PREQUAL for each set orthogroup is  easy:

```
for f in *fa; do prequal $f ; done
```

The filtered (masked) alignments are in `.filtered` whereas `.prequal` contains relevant information such as the number of residues filtered.



## Multiple sequence alignment


The next step is to infer multiple sequence alignments. Multiple sequence alignments allow us to *know* which amino acids/ nucleotides are homologous. A simple yet accurate tool is [MAFFT](https://mafft.cbrc.jp/alignment/server/).

We will align gene files separately using a for loop:

```
for f in *filtered; do mafft $f > $f.mafft; done
```


## Alignment trimming


Some gene regions (e.g., fast-evolving) are difficult to align and thus positional homology is uncertain. It is unclear (probably problem-specific) whether trimming badly-aligned regions [improves](https://academic.oup.com/sysbio/article/56/4/564/1682121) or [worsens](https://academic.oup.com/sysbio/article/64/5/778/1685763) tree inference. However, gently trimming very incomplete positions (e.g. with >80% gaps) will speeds up computation in the next steps without significant information loss.

To trim alignment positions we can use [BMGE](https://bmcevolbiol.biomedcentral.com/articles/10.1186/1471-2148-10-210) but several other software are also available.

To remove alignment positions with > 80% gaps:

```
for f in *mafft; do java -jar /software/BMGE.jar -i $f -t AA -g 0.8 -h 1 -w 1 -of $f.g08.fas; done
```

Alternatively, the default settings in BMGE will remove incomplete positions and additionally trim high-entropy (likely fast-evolving) positions:

```
for f in *mafft; do java -jar /software/BMGE.jar -i $f -t AA -of $f.bmge.fas; done
```

While diving into phylogenomic pipelines, it is always advisable to check a few intermediate results to ensure we are doing what we should be doing. Multiple sequence alignments can be visualized in [SeaView](http://doua.prabi.fr/software/seaview) or [AliView](https://github.com/AliView/AliView). Also, one could have a quick look at alignments using command line tools (`less -S`).


## Concatenate alignment


To infer our phylogenomic tree we need to concatenate single-gene alignments. This can be done with tools such as [FASconCAT](https://github.com/PatrickKueck/FASconCAT-G), which will read in all `\*.fas` `\*.phy` or `\*.nex` files in the working directory and concatenate them (in random order).

```
perl /software/FASconCAT-G_v1.04.pl -l -s
```

Is your concatenated file what you expected? It should contain 23 taxa and 21 genes. You might check the concatenation (`FcC_supermatrix.fas`) and the file containing the coordinates for gene boundaries (`FcC_supermatrix_partition.txt`). Looking good? Then your concatenated dataset is ready to rock!!



## Concatenation: Maximum likelihood


One of the most common approaches in phylogenomics is gene concatenation: the signal from multiple genes is "pooled" together with the aim of increasing resolution power. This method is best when among-gene discordance is low.

We will use [IQTREE](http://www.iqtree.org/), an efficient and accurate software for maximum likelihood analysis. Another great alternative is [RAxML](https://github.com/stamatak/standard-RAxML). The most simple analysis is to treat the concatenated dataset as a single homogeneous entity. We need to provide the number of threads to use (`-nt 1`) input alignment (`-s`), tell IQTREE to select the best-fit evolutionary model with BIC (`-m TEST -merit BIC -msub nuclear`) and ask for branch support measures such as non-parametric bootstrapping and approximate likelihood ratio test (`-bb 1000 -alrt 1000 -bnni`):

```
iqtree -s FcC_supermatrix.fas -m TEST -msub nuclear -bb 1000 -alrt 1000 -nt 1 -bnni -pre unpartitioned
```

A more sophisticated approach would be to perform a partitioned maximum likelihood analysis, where different genes (or other data partitions) are allowed to have different evolutionary models. This should provide a better fit to the data but will increase the number of parameters too. To launch this analysis we need to provide a file containing the coordinates of the partitions (`-spp`) and we can ask IQTREE to select the best-fit models for each partition, in this case according to AICc (more suitable for shorter alignments).

```
iqtree -s FcC_supermatrix.fas -spp FcC_supermatrix_partition.txt -m TEST -msub nuclear -merit AICc -bb 1000 -alrt 1000 -nt 1 -bnni -pre partitioned
```

Congratulations!! If everything went well, you should get your maximum likelihood estimation of the vertebrate phylogeny (`.treefile`)! Looking into the file you will see a tree in parenthetical (newick) format. See below how to create a graphical representation of your tree.



## Coalescence analysis


An alternative to concatenation is to use a multispecies coalescent approach. Unlike maximum likelihood, coalescent methods account for incomplete lineage sorting (ILS; an expected outcome of evolving populations). These methods are particularly useful when we expect high levels of ILS, e.g. when speciation events are rapid and leave little time for allele coalescence.

We will use [ASTRAL](https://github.com/smirarab/ASTRAL), a widely used tool that scales up well to phylogenomic datasets. It takes a set of gene trees as input and will generate the coalescent "species tree". ASTRAL assumes that gene trees are estimated without error.

Thus, before running ASTRAL, we will need to estimate individual gene trees. This can be easily done calling IQTREE in a for loop:

```
for f in *filtered.mafft.g08.fas; do iqtree -s $f -m TEST -msub nuclear -merit AICc -nt 1; done
```

After all gene trees are inferred, we should put them all into a single file:

```
cat *filtered.mafft.g08.fas.treefile > my_gene_trees.tre
```

Now running ASTRAL is trivial, providing the input file with the gene trees and the desired output file name:

```
java -jar /srv/evop/resources/astral/ASTRAL/Astral/astral.5.7.3.jar -i my_gene_trees.tre -o species_tree_ASTRAL.tre 2> out.log
```

Congratulations!! You just got your coalescent species tree!! Is it different from the concatenated maximum likelihood trees? 



## Tree visualization


Trees are just text files representing relationships with parentheses; did you see that already? But it is more practical to plot them as a graph, for which we can use tools such as [iTOL](https://itol.embl.de), [FigTree](https://github.com/rambaut/figtree/releases), [iroki](https://www.iroki.net/), or R (e.g. [ggtree](https://bioconductor.org/packages/release/bioc/html/ggtree.html), [phytools](http://www.phytools.org/); see provided R code).

Upload your trees to iTOL. Trees need to be rooted with an outgroup. Click in the branch of *Callorhinchus milii* and the select "Tree Structure/Reroot the tree here". Branch support values can be shown under the "Advanced" menu. The tree can be modified in many other ways, and finally, a graphical tree can be exported. Similar options are available in FigTree.

[Well done!](https://media.giphy.com/media/wux5AMYo8zHgc/giphy.gif)

