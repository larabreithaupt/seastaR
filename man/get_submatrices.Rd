\name{get_submatrices}
\alias{get_submatrices}

\title{
get_submatrices
}
\description{
    This function builds the full variance-covariance matrix by using theory to calculate submatrices for each triplet.
}
\usage{
    get_submatrices(sptree)
}

\arguments{
  \item{sptree}{
     Species tree phylo object produced from file parser (parse_input_file); species tree is in coalescent units
  }
}
\value{
  returns a list of all of the submatrices of the species tree
}
\references{
  NEED TO ADD
}
\author{
    Lara Breithaupt and Mark Hibbins
}

\examples{
    sub_matrices <- get_submatrices(test_input_files/seastar_sptree_test_input.txt)
}
