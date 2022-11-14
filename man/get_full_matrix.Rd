\name{get_full_matrix}
\alias{get_full_matrix}

\title{
get_full_matrix
}
\description{
    This function builds a covariance matrix from an input tree using expectations from the multispecies coalescent
}
\usage{
    get_full_matrix(sptree)
}

\arguments{
  \item{sptree}{
     Species tree phylo object. Must be ultrametric and in coalescent units
    }
}
\value{
  Returns a ntaxa x ntaxa matrix containing variance and covariance values for each pair of species
}
\references{
  Revell LJ (2012). “phytools: An R package for phylogenetic comparative biology (and other things).” Methods in Ecology and Evolution, 3, 217-223.
}
\author{
    Lara Breithaupt and Mark Hibbins
}

}
\examples{
  sptree <- parse_input_file(test_input_files/seastar_sptree_test_input.txt, genetrees = FALSE)
  full_matrix <- get_full_matrix(sptree)
}
