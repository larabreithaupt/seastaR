remove(list=ls())
library(phytools)
library(geiger)
library(tidyverse)
lib <- modules::use("R")

#Load in files with my paths (Lara - comment these out on your machine)
#tree_list <- lib$tools$parse_input_file("C:/Users/18126/OneDrive - Indiana University/Projects/pruning_alg/seastaR/test_input_files/seastar_genetrees_test_input.txt", genetrees = TRUE)
#sptree <- lib$tools$parse_input_file("C:/Users/18126/OneDrive - Indiana University/Projects/pruning_alg/seastaR/test_input_files/seastar_sptree_test_input.txt", genetrees = FALSE)

#Load in files for genetrees and sptree
tree_list <- lib$tools$parse_input_file("~/Downloads/Hahn Lab Pruning Algorithm/seastar/seastaR/test_input_files/seastar_genetrees_test_input.txt", genetrees = TRUE)
#sptree <- lib$tools$parse_input_file("~/Downloads/Hahn Lab Pruning Algorithm/seastar/seastaR/test_input_files/seastar_sptree_test_input.txt", genetrees = FALSE)
#sptree <- lib$tools$parse_input_file("~/Downloads/Hahn Lab Pruning Algorithm/seastar/seastaR/test_input_files/three_taxa_tree_test.txt", genetrees = FALSE)
sptree <- lib$tools$parse_input_file("~/Downloads/Hahn Lab Pruning Algorithm/seastar/seastaR/test_input_files/symm_four_taxa_tree_test.txt", genetrees = FALSE)

#Get genetree matrix
genetree_vcv <- lib$estimated_trees$trees_to_vcv(tree_list)
#print(genetree_vcv)

#Get sptree matrix
full_matrix <- lib$sptree_theory$get_full_matrix(sptree)
#print(full_matrix)

#Get sigma squared and simulate traits

trait <- lib$inference$simulate_traits(1000, full_matrix, 1)
sigma_2_est <- apply(trait, 1, lib$inference$sigma2_inference, var_covar = full_matrix)
sigma_2 <- lib$inference$sigma2_inference(full_matrix, trait)



