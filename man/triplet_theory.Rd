\name{triplet_theory}
\alias{triplet_theory}

\title{
triplet_theory
}
\description{
    This function builds a 3x3 submatrix using the multispecies coalescent on the input triplet
}
\usage{
    triplet_theory(tau, height, triplet)
}

\arguments{
  \item{tau}{
     The triplet internal branch in coalescent units
  }
  \item{height}{
     The height of the triplet MRCA in coalescent units
  }
  \item{triplet}{
     A character vector of length 3, containing the name of each taxon in the triplet. Sister taxa must be listed first.
  }
}
\details{
    When using get_full_matrix, this function is called internally to construct the partial covariance matrix for each triplet.
}
\value{
    Returns a partial matrix that will be used to construct the final matrix for all species
}
\references{
    Revell LJ (2012). “phytools: An R package for phylogenetic comparative biology (and other things).” Methods in Ecology and Evolution, 3, 217-223.
}
\author{
    Lara Breithaupt and Mark Hibbins
}

\examples{
    sub_matrix <- triplet_theory(0.87, 2.2, c("sp4", "sp3", "sp2"))
}
