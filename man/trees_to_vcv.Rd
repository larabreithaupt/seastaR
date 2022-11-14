\name{trees_to_vcv}
\alias{trees_to_vcv}

\title{
trees_to_vcv
}
\description{
    Builds a phylogenetic variance/covariance matrix from a list of input gene trees. 
}
\usage{
    trees_to_vcv(tree_list)
}

\arguments{
  \item{tree_list}{
      List of phylo objects; last element of list must be a vector of genetree frequencies. Gene trees must be ultrametric. 
  }
}
\details{
    This function will take a list of genetrees and call get_partial_cov_matrix to get the partial matrices from each genetree to construct the final matrix.
}
\value{
    Returns an n x n matrix where n is the number of taxa.
}
\references{
    Revell LJ (2012). “phytools: An R package for phylogenetic comparative biology (and other things).” Methods in Ecology and Evolution, 3, 217-223.
}
\author{
    Lara Breithaupt and Mark Hibbins
}

\examples{
  genetree_list <- parse_input_file(test_input_files/seastar_genetrees_test_input.txt)
  genetree_cov_matrix <- trees_to_vcv(genetree_list)
}

