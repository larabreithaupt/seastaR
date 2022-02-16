#import(utils)

#Functions for building covariance matrix from input set of
#estimated gene trees

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
