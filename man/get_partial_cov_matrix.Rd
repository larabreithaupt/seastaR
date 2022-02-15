\name{get_partial_cov_matrix}
\alias{get_partial_cov_matrix}

\title{
get_partial_cov_matrix
}
\description{
    Gets a partially filled matrix with covariances from one gene tree
}
\usage{
    get_partial_cov_matrix(genetree)
}

\arguments{
  \item{genetree}{
      A gene tree phylo object. Must be ultrametric. 
  }
}
\details{
    When trees_to_vcv is used, this function is called internally to fill the part of the matrix corresponding to each gene tree in the input list.
}
\value{
    Returns a phylogenetic variance-covariance matrix of size #taxa x #taxa
}
\references{
    Revell LJ (2012). “phytools: An R package for phylogenetic comparative biology (and other things).” Methods in Ecology and Evolution, 3, 217-223.
}
\author{
    Lara Breithaupt and Mark Hibbins
}
\examples{
    genetree_list <- parse_input_file(test_input_files/seastar_genetrees_test_input.txt)
    partial_matrix <- get_partial_cov_matrix(genetree_list[[1]])
}

