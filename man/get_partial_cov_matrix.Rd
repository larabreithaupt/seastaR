\name{get_partial_cov_matrix}
\alias{get_partial_cov_matrix}

\title{
get_partial_cov_matrix
}
\description{
    Gets a partially filled matrix with covariances from one genetree
}
\usage{
    get_partial_cov_matrix(genetree)
}

\arguments{
  \item{genetree}{
      Singular genetree is a phylogenetic object
  }
}
\details{
    this function will be called on each genetree in the list of genetree phylogenetic objects to construct the full matrix
}
\value{
    returns a phylogenetic variance-covariance or correlation matrix of size #taxa x #taxa
}
\references{
    NEED TO ADD
}
\author{
    Lara Breithaupt and Mark Hibbins
}
\examples{
    partial_matrix <- get_partial_cov_matrix(test_input_files/seastar_genetrees_test_input.txt[[i]])
}

