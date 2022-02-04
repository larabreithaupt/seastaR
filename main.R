remove(list=ls())
library(phytools)
library(geiger)
library(tidyverse)
lib <- modules::use("R")

#Load in files with my paths (Lara - comment these out on your machine)
#tree_list <- lib$tools$parse_input_file("C:/Users/18126/OneDrive - Indiana University/Projects/pruning_alg/seastaR/test_input_files/seastar_genetrees_test_input.txt", genetrees = TRUE)
#sptree <- lib$tools$parse_input_file("C:/Users/18126/OneDrive - Indiana University/Projects/pruning_alg/seastaR/test_input_files/seastar_sptree_test_input.txt", genetrees = FALSE)

tree_list <- lib$tools$parse_input_file("~/Downloads/Hahn Lab Pruning Algorithm/seastar/seastaR/test_input_files/seastar_genetrees_test_input.txt", genetrees = TRUE)
sptree <- lib$tools$parse_input_file("~/Downloads/Hahn Lab Pruning Algorithm/seastar/seastaR/test_input_files/seastar_sptree_test_input.txt", genetrees = FALSE)


triplet_branch <- lib$sptree_theory$get_internal_branch(c("sp1", "sp2", "sp3"), sptree)

triplets <- lib$sptree_theory$get_triplet_branches(sptree)
final_matrix <- lib$sptree_theory$get_submatrices(sptree)

#This works fine
genetree_vcv <- lib$estimated_trees$trees_to_vcv(tree_list)

full_matrix <- lib$sptree_theory$get_full_matrix(sptree)


