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
  
  #in-order tree traversal
  
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
  
  tree1 = tree_list[[1]]
  tips = tree1[["tip.label"]]
  len_tip = length(tips)
  
  genetree_vcv <- matrix(0, 
                        len_tip, 
                        len_tip)
  rownames(genetree_vcv) <- tips
  colnames(genetree_vcv) <- tips

  for(i in 1:(length(tree_list) - 1)){
    
    tree <- tree_list[[i]]
  
    edge_len <- length(tree[["edge.length"]])
    
    height <- max(nodeHeights(tree))
    
    tip_combos <- combn(tips, 2)
    
    for(j in 1:length(tip_combos[1,])){
      
      col <- tip_combos[,j]
      
      mrca <- findMRCA(tree, col)
      
      if (mrca == (len_tip + 1)){
        next
      }
      
      else{
        
        tree_edge <- tree[["edge"]]
        tree_edge_len <- tree[["edge.length"]]
        
        for(k in 1:length(tree_edge_len)){
          
          if((tree_edge[k, 1] == (mrca - 1)) && (tree_edge[k, 2] == mrca)){
            
            internal <- tree_edge_len[k]
            gene_freq <- tree_list[[5]][i]
            
            add1 <- (genetree_vcv[col[1], col[2]]) + (internal*gene_freq)
            add2 <- (genetree_vcv[col[2], col[1]]) + (internal*gene_freq)
            
            genetree_vcv[col[1], col[2]] <- add1
            genetree_vcv[col[2], col[1]] <- add2
            
          }
          
        }
        
        for(m in 1:len_tip){
          
          label <- tips[m]
          diag <- genetree_vcv[label, label] + (height*gene_freq)
          
          genetree_vcv[label, label] <- diag
          
        }
          
      }
      
    }  
    
  }
  print(genetree_vcv)
  return(genetree_vcv)
  
}

genetree_vcv <- trees_to_vcv(tree_list)

