\name{get_triplet_branches}
\alias{get_triplet_branches}

\title{
get_triplet_branches
}
\description{
    This function returns a data frame containing the internal branch lengths and root node heights for all triplets contained in the species tree
}
\usage{
    get_triplet_branches(sptree)
}

\arguments{
  \item{sptree}{
     Species tree phylo object. Must be ultrametric and in coalescent units. 
  }
}
\value{
    Returns a data frame with the triplet name, its internal branch length, and the height of its MRCA
}
\references{
    Revell LJ (2012). “phytools: An R package for phylogenetic comparative biology (and other things).” Methods in Ecology and Evolution, 3, 217-223.
}
\author{
    Lara Breithaupt and Mark Hibbins
}

\examples{
    sptree <- parse_input_file(test_input_files/seastar_sptree_test_input.txt, genetrees = FALSE)
    branches_df <- get_triplet_branches(sptree)
}

