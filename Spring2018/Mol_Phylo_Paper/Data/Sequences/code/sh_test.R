setwd("D:/GitHub/graduateschoolclasses/Spring2018/Mol_Phylo_Paper/Data/Sequences/trials/final/")

library(phangorn)
library(phylobase)
library(phytools)

#Conversion of NEXUS format to phyDat (Matrix)
three_gene_seqs <- read.phyDat(file = "three_gene_aligned.fas", format = "fasta")
  two_gene_seqs <- read.phyDat(file = "Lu_et_al_2014.fst", format = "fasta")

#Assigning trees to a variable
three_gene_tree <-  read.newick(file = "three_gene_tree.newick")
  two_gene_tree <- read.newick(data = "Lu_et_al_2014.newick")

#Likelihood computing for a tree
#Have to run in order to get likelihood into R
three_gene_tree_fit <- pml(three_gene_tree, three_gene_seqs)
  two_gene_tree_fit <- pml(two_gene_tree, two_gene_seqs)

#Shimodaira-Hasegawa test to see if trees are significantly different  
SH_test_for_trees <- SH.test(three_gene_tree_fit, two_gene_tree_fit, B = 100000)

SH_test_for_trees






