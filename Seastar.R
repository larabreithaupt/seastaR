remove(list=ls())
library(ape)
library(phytools)
library(geiger)
library(tidyverse)

parse_input_file <- function(file_path, genetrees = TRUE){

  ### Reads input file, containing 
  ### gene trees in Newick format and
  ### their frequencues
  

  
  if (genetrees){
  
    iftree = FALSE
    tree_list <- vector(mode = "list")
    freqs_vector <- vector()
    input_file = file(file_path, "r")
    
    i <- 1
    j <- 1
  
    while (TRUE) {
    
      line = readLines(input_file, n = 1)
      if (length(line) == 0){ #if no more lines, break out of loop
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
    
  } else{
      
    input_file = file(file_path, "r")
    
    text_line <- readLines(input_file, n = 1)
    tree <- read.tree(text = text_line)

    close(input_file)
    return(tree)
  }
  
}

tree_list <- parse_input_file("~/Downloads/Hahn Lab Pruning Algorithm/seastar/seastaR/test_input_files/seastar_genetrees_test_input.txt", genetrees = TRUE)
sptree <- parse_input_file("~/Downloads/Hahn Lab Pruning Algorithm/seastar/seastaR/test_input_files/seastar_sptree_test_input.txt", genetrees = FALSE)

get_internal_branch <- function(triplet, sptree){
  
  triplet_mrca <- phytools::findMRCA(sptree, triplet)
  print(triplet_mrca)
  tips_pairwise <- combn(triplet, 2)
  print(tips_pairwise)
  
  sister <- ""
  sister_mrca <- 0
    
  for(i in 1:length(tips_pairwise[1,])){#for each pairwise combo
    
    combo <- tips_pairwise[,i]
    
    mrca <- phytools::findMRCA(sptree, combo)
    
    if(mrca > triplet_mrca){
      
      sister <- combo
      sister_mrca <- mrca
      
    }
    
  }
  print(sister_mrca)
  
}

internal <- get_internal_branch(c("sp4", "sp3", "sp1"), sptree)

get_quartets <- function(sptree){
  #pulls out the quartets and internal branch lengths
  
  
  
}

quartet <- get_quartets(tree_list)


get_theory_matrix <- function(sptree){
  
  
  
}

get_partial_cov_matrix <- function(genetree){
  
  #Gets a partially filled matrix with covariances from one tree 
  
  tree_height <- max(phytools::nodeHeights(genetree))
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



