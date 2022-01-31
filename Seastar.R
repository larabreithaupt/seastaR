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
  tips_pairwise <- combn(triplet, 2)
  
  internal <- 0

  
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
  tree_edge <- sptree[["edge"]]
  tree_edge_len <- sptree[["edge.length"]]
  
  counter <- triplet_mrca
  
  for(j in 1:length(tree_edge_len)){ #iterating through node traversal
    
    
    if (counter == sister_mrca){
      
      break
      
    }
      
    else if ((tree_edge[j, 1] == (counter)) && (tree_edge[j, 2] == (counter + 1))){ #traversal from parent node to our MRCA
      
      internal <- internal + tree_edge_len[j] #internal branch length
      
      counter <- counter + 1
    }
  }
  
  tree_height <- max(phytools::nodeHeights(sptree))
  
  return(list(internal, tree_height, sister))
}

internal <- get_internal_branch(c("sp1", "sp2", "sp3"), sptree)

get_triplet_branches <- function(sptree){
  #pulls out the quartets and internal branch lengths
  
  triplet_internals <- list()
  
  tip_labels <- sptree[["tip.label"]]
  triplets <- combn(tip_labels, 3)
  
  for(i in 1:length(triplets[1, ])){
    
    internal_branch <- get_internal_branch(triplets[, i], sptree)
    
    len <- internal_branch[1]
    sisters <- internal_branch[2]
    triplet <- sisters[[1]]
    
    print(sisters)
    
    for (j in 1:length(triplets[, i])){
      
      if (!(triplets[j, i] %in% triplet)) {
        triplet[3] <- triplets[j, i] 
      }
      
    }
    
    tip_concat <- paste(triplet, collapse = "")
    
    triplet_internals[tip_concat] = len
    
  }

  return(triplet_internals)

}

triplets <- get_triplet_branches(sptree)


triplet_theory <- function(tau, height, triplet) {
  
  partial_matrix <- matrix(c(0, 0, 0,
                             0, 0, 0, 
                             0, 0, 0), nrow = 3, ncol = 3)
 
 # print(triplet)
  rownames(partial_matrix) <- triplet
  colnames(partial_matrix) <- triplet
  
  LS <- 1 - exp(-tau)
  ILS <- (1/3)*exp(-tau)
  
  AB_covar <- LS*(tau + (tau/(exp(tau)-1))) + ILS
  
  #Discordant trees have internal branch length of 1
  BC_covar <- ILS
  AC_covar <- ILS
  
  #Weighted height of all gene trees 
  all_var <- LS*(height + 1) + 3*ILS*(height + 1 + 1/3)
  
  
  #Fill matrix 
  
  partial_matrix[1,1] <- all_var
  partial_matrix[2,2] <- all_var
  partial_matrix[3,3] <- all_var
  
  partial_matrix[1,2] <- AB_covar
  partial_matrix[2,1] <- AB_covar
  
  partial_matrix[1,3] <- AC_covar
  partial_matrix[3,1] <- AC_covar
  
  partial_matrix[2,3] <- BC_covar
  partial_matrix[3,2] <- BC_covar
  
  return(partial_matrix)
}


get_theory_matrix <- function(sptree){
  
  tips = sptree[["tip.label"]]
  len_tip = length(tips)
  
  branches <- get_triplet_branches(sptree)
  
  #Initialize Matrix
  theory_matrix <- matrix(0, len_tip, len_tip)
  rownames(theory_matrix) <- tips
  colnames(theory_matrix) <- tips
  
  for (i in 1:length(branches[1])){
    
    sub_matrix <- triplet_theory(branches[1], branches[2], branches[3])
    
    #print(sub_matrix)
  }
  
}

final_matrix <- get_theory_matrix(sptree)

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
#print(genetree_vcv)



