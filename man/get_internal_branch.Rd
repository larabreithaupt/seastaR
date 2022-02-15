\name{get_internal_branch}
\alias{get_internal_branch}

\title{
get_internal_branch
}
\description{
    Function to get the internal branch length contained in a single specified
    triplet, which is pulled out of the larger species tree
}
\usage{
    get_internal_branch(triplet, sptree)
}

\arguments{
  \item{triplet}{
    A character vector of length 3 containing the name of each taxon in the triplet. The sister taxa in the tripet should be listed first.
}
  \item{sptree}{
    Species tree phylo object. Must be ultrametric and in coalescent units.
  }
}
\details{
    When get_full_matrix is used, this function is called internally on all triplets contained in the species tree to get their internal branches. 
}
\value{
    A list with three elements: the length of the internal branch, the node height of the triplet MRCA, and a character vector of length 2 containing the names of the sister taxa in the triplet.
}
\references{
    Revell LJ (2012). “phytools: An R package for phylogenetic comparative biology (and other things).” Methods in Ecology and Evolution, 3, 217-223.
}
\author{
    Lara Breithaupt and Mark Hibbins
}
\examples{
  sptree <- parse_input_file(test_input_files/seastar_sptree_test_input.txt, genetrees = FALSE)
  internal_branch <- get_internal_branch(c("sp3", "sp4", "sp1"), sptree)
}

