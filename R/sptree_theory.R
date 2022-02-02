import(tidyverse)
import(utils)
import(stringr)

#Functions for building covariance matrix from input 
#species tree in coalescent units

get_internal_branch <- function(triplet, sptree){
  
  #This function gets the internal branch length contained 
  #in a single specified triplet, which is pulled out of the
  #larger species tree
  
  triplet_mrca <- phytools::findMRCA(sptree, triplet)
  
  max_height <- max(phytools::nodeHeights(sptree))
  triplet_height <- phytools::nodeheight(sptree, triplet_mrca)
  final_height <- max_height - triplet_height
  
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
  
  return(list(internal, final_height, sister))
}

get_triplet_branches <- function(sptree){
  
  #This function returns a list of internal branch lengths
  #for all triplets contained in the species tree
  
  tau <- c()
  height <- c()
  triplet_names <- c()
  
  tip_labels <- sptree[["tip.label"]]
  triplets <- combn(tip_labels, 3)
  
  for(i in 1:length(triplets[1, ])){
    
    internal_branch <- get_internal_branch(triplets[, i], sptree)
    #print(internal_branch)
    
    len <- internal_branch[1]
    sisters <- internal_branch[3]
    triplet <- sisters[[1]]
    
    for (j in 1:length(triplets[, i])){
      
      if (!(triplets[j, i] %in% triplet)) {
        triplet[3] <- triplets[j, i] 
      }
      
    }
    
    tip_concat <- paste(triplet, collapse = "_")
    
    tau <- append(tau, len)
    triplet_names <- append(triplet_names, tip_concat)
    height <- append(height, internal_branch[2])
    
  }
  
  ret <- as.data.frame(cbind(tau, height, triplet_names))
  return(ret)
}

triplet_theory <- function(tau, height, triplet) {
  
  #This function builds a 3x3 submatrix constructed 
  #using the multispecies coalescent on the input 
  #triplet
  
  partial_matrix <- matrix(c(0, 0, 0,
                             0, 0, 0, 
                             0, 0, 0), nrow = 3, ncol = 3)
  
  # print(triplet)
  names <- strsplit(as.character(triplet), split = "_")[[1]]
  rownames(partial_matrix) <- names
  colnames(partial_matrix) <- names
  
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
  
  #This function builds the full variance-covariance matrix
  #by using theory to calculate submatrices for each triplet.
  #In progress
  
  tips = sptree[["tip.label"]]
  len_tip = length(tips)
  
  branches_df <- get_triplet_branches(sptree)
  #print(branches_df)
  
  #Initialize Matrix
  theory_matrix <- matrix(0, len_tip, len_tip)
  rownames(theory_matrix) <- tips
  colnames(theory_matrix) <- tips
  
  for (i in 1:length(branches_df[1, ])){
    
    sub_matrix <- triplet_theory(as.numeric(branches_df[i, 1]), 
                                 as.numeric(branches_df[i, 2]), 
                                 branches_df[i, 3])
    
    print(sub_matrix)
  }
  
}