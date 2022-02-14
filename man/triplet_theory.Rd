\name{triplet_theory}
\alias{triplet_theory}

\title{
triplet_theory
}
\description{
    This function builds a 3x3 submatrix constructed using the multispecies coalescent on the input triplet
}
\usage{
    triplet_theory(tau, height, triplet)
}

\arguments{
  \item{tau}{
     the divergence time measured in expected number of mutations of the triplet; used to compute Lineage Sorting and Incomplete          Lineage Sorting to calculate variances and covariances
  }
  \item{height}{
     height of the triplet; used to compute Lineage Sorting and Incomplete Lineage Sorting to calculate variances and covariances
  }
  \item{triplet}{
     triplet from the species tree
  }
}
\details{
    the input values will be pulled from the data frame produced in get_triplet_branches, and called on each value in the dataframe      that each include a tau value, a height, and a triplet
}
\value{
    returns a partial matrix that will be used to construct the final species tree matrix
}
\references{
    NEED TO ADD
}
\author{
    Lara Breithaupt and Mark Hibbins
}

\examples{
    sub_matrix <- triplet_theory(as.numeric(branches_df[i, 1]),
                                 as.numeric(branches_df[i, 2]),
                                 branches_df[i, 3])
}
