remove(list=ls())
library(ape)
library(phytools)
library(geiger)

parse_input_file <- function(file_path){

  ### Reads input file, containing 
  ### gene trees in Newick format and
  ### their frequencues
  
  iftree =FALSE
  tree_list <- vector(mode = "list")
  freqs_vector <- vector()
  input_file = file(file_path, "r")
  
  i <- 1
  j <- 1
  
  while (TRUE) {
    
    line = readLines(input_file, n = 1)
    if (length(line) == 0){ #skip empty lines
      break
    }
    
    if (startsWith(line, "(")){ #line containing a newick tree
      tree <- read.tree(text = line)
      tree_list[[i]] <- tree
      i <- i + 1
    }
    
    if (startsWith(line, "0")){ #line containing tree frequency 
      freq <- as.double(line)
      freqs_vector[j] <- freq
      j <- j + 1
    }
    
    tree_list[[i]] <- freqs_vector #add vector of tree freqs to end of list 
    
  }

  close(input_file)
  
  return(tree_list)
  
}

tree_list <- parse_input_file("~/Downloads/Hahn Lab Pruning Algorithm/seastar/seastaR/seastar_test_input.txt")

trees_to_phylo_star <- function(genetrees){

  #Function that takes list of gene trees as input 
  #and returns a vcv and single phylo object summarizing them.
  #Currently hard-codes an example case, need to create
  #a generalized version. This function is the main function
  #of the program. 
  
  combined_trees <- list()
  class(combined_trees) <- "phylo"
  
  combined_trees$edge <- matrix(c(
    4,5, #Pre-order tree traversal that contains extra
    5,6, #nodes for covariances not described by a standard species 
    6,3, #tree
    6,1,
    5,7,
    7,2,
    7,1, 
    4,1), 8,2, byrow = TRUE)
  combined_trees$Nnode <- 4 
  combined_trees$node.label <- c(4,5,6,7) #internal node labels 
  combined_trees$tip.label <- c("A", "B", "C") #leaf node labels 
  combined_trees$edge.length <- c(0.5, 0.1, 0.4, 0.4, 0.3, 0.2, 0.2, 1) #branch lengths corresponding to traversals 
  
  combined_vcv <- ape::vcv(combined_trees) #variance-covariance matrix from our phylo object
  
  return(list(combined_vcv, combined_trees))
}

#new_phylo <- trees_to_phylo_star(tree_list)
#print(new_phylo[[1]])
#print(det(new_phylo[[1]]))


#tips <- c(1, 2, 3)
#names(tips) <- c("A", "B", "C")

#anc <- geiger::fitContinuous(new_phylo[[2]], tips)
#print(anc)

get_partial_cov_matrix <- function(genetree){
  
  #Gets a partially filled matrix with covariances from one tree 
  
  tree_height <-max(phytools::nodeHeights(genetree))
  tips = genetree[["tip.label"]]
  len_tip = length(tips)
  
  tip_combos <- combn(tips, 2) #pairwise combinations of taxa
  
  #Initialize matrix
  genetree_vcv <- matrix(0, 
                         len_tip, 
                         len_tip)
  rownames(genetree_vcv) <- tips
  colnames(genetree_vcv) <- tips
  
  for(j in 1:length(tip_combos[1,])){ #for each pairwise combo
    
    col <- tip_combos[,j] #current pairwise combo
    mrca <- phytools::fastMRCA(genetree, col[1], col[2]) #most recent common ancestor of combo
    mrca_height <- phytools::fastHeight(genetree, col[1], col[2]) #height of mrca from root
    
    genetree_vcv[col[1], col[2]] <- mrca_height
    genetree_vcv[col[2], col[1]] <- mrca_height
    
  }
  
  for(m in 1:len_tip){ #Fill the diagonal elements (variance of each species)
    
    label <- tips[m]
    genetree_vcv[label, label] <- tree_height
    
  }
  
  return(genetree_vcv)
  
}

trees_to_vcv <- function(tree_list) {

  #Our own implementation for constructing a 
  #phylogenetic covariance matrix from a set of 
  #input gene trees. Used to benchmark trees_to_star,
  #and also for methods where you can pass the covariance
  #matrix directly. 
  
  tree1 = tree_list[[1]]
  tips = tree1[["tip.label"]] 
  len_tip = length(tips)
  
  #Initialize vcv
  genetree_vcv <- matrix(0, 
                        len_tip, 
                        len_tip)
  rownames(genetree_vcv) <- tips
  colnames(genetree_vcv) <- tips
  

  for(i in 1:(length(tree_list) - 1)){
    
    gene_freq <- tree_list[[length(tree_list)]][i] #gene tree frequency 
    
    partial_matrix <- get_partial_cov_matrix(tree_list[[i]]) #get partial matrix
    
    scaled_matrix <- partial_matrix*gene_freq #scale covs
    scaled_matrix <- scaled_matrix[, order((colnames(scaled_matrix)))] #reorder cols
    scaled_matrix <- scaled_matrix[order((rownames(scaled_matrix))), ] #reorder ros
    
    genetree_vcv <- genetree_vcv + scaled_matrix
    
  }
  
  return(genetree_vcv)  
}

genetree_vcv <- trees_to_vcv(tree_list)
print(genetree_vcv)



