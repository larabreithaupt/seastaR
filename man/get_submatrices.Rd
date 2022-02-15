\name{get_submatrices}
\alias{get_submatrices}

\title{
get_submatrices
}
\description{
    This function builds a list of 3x3 variance-covariance submatrices from all triplets contained in the input species tree.
}
\usage{
    get_submatrices(sptree)
}

\arguments{
  \item{sptree}{
     Species tree phylo object. Must be ultrametric and in coalescent units. 
  }
}
\value{
  Returns a list of all of the submatrices of the species tree
}
\references{
  Revell LJ (2012). “phytools: An R package for phylogenetic comparative biology (and other things).” Methods in Ecology and Evolution, 3, 217-223.
}
\author{
    Lara Breithaupt and Mark Hibbins
}

\examples{
    sptree <- parse_input_file(test_input_files/seastar_sptree_test_input.txt, genetrees = FALSE)
    sub_matrices <- get_submatrices(sptree)
}
