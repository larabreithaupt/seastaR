remove(list=ls())
library(ape)
library(phytools)
library(geiger)

parse_input_file <- function(file_path){
  
  iftree = FALSE
  tree_list <- vector(mode = "list")
  freqs_vector <- vector()
  input_file = file(file_path, "r")
  
  i <- 1
  j <- 1
  
  while (TRUE) {
    
    line = readLines(input_file, n = 1)
    if (length(line) == 0){
      break
    }
    
    if (startsWith(line, "(")){
      tree <- read.tree(text = line)
      tree_list[[i]] <- tree
      i <- i + 1
    }
    
    if (startsWith(line, "0")){
      freq <- as.double(line)
      freqs_vector[j] <- freq
      j <- j + 1
    }
    
    tree_list[[i]] <- freqs_vector
    
  }

  close(input_file)
  
  return(tree_list)
  
}

tree_list <- parse_input_file("~/Downloads/Hahn Lab Pruning Algorithm/seastar/seastar_test_input.txt")

trees_to_star <- function(genetrees){
  
  combined_trees <- list()
  class(combined_trees) <- "phylo"
  
  combined_trees$edge <- matrix(c(
    4,5,
    5,6,
    6,3,
    6,1,
    5,7,
    7,2,
    7,1, 
    4,1), 8,2, byrow = TRUE)
  combined_trees$Nnode <- 4
  combined_trees$node.label <- c(4,5,6,7)
  combined_trees$tip.label <- c("A", "B", "C")
  combined_trees$edge.length <- c(0.5, 0.1, 0.4, 0.4, 0.3, 0.2, 0.2, 1)
  
  #trait <- setNames(c(0.2, 0.2, 0.2, 0.2), c("A", "B", "C"))
  #print(trait)
  #fit.BM <- ace(trait, combined_trees, method = "ML")
  
  tree <- vcv(combined_trees)
  
  return(list(tree, combined_trees))
  
}

new_phylo <- trees_to_star(tree_list)

trees_to_vcv <- function(tree_list) {
  
  len_tip = length(tree_list[[1]]$tip.label)
  
  genetree_vcv <- matrix(0, 
                        len_tip, 
                        len_tip)
  
  for(i in 1:length(tree_list)){
    
  
    tree <- tree_list[[i]]
    
    edge_len <- length(tree$edge.length)
    
    
    for(j in 1:edge_len){
      
      node1 = tree$edge[j, 1]
      node2 = tree$edge[j, 2]
      
      if ((node1 > len_tip) && (node2 > len_tip)) {
        
        ind = tree$edge.length[j]
        
        print("internal branch")
        
      }
      
    }
    
  }
  
  return(genetree_vcv)
  
}

genetree_vcv <- trees_to_vcv(tree_list)

