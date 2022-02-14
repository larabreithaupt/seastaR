\name{trees_to_vcv}
\alias{trees_to_vcv}

\title{
trees_to_vcv
}
\description{
    Our implementation for constructing a phylogenetic covariance matrix from a set of input gene trees
}
\usage{
    trees_to_vcv(tree_list)
}

\arguments{
  \item{tree_list}{
      List of gene tree phylogenetic objects produced from file parser (parse_input_file); last element of list is a vector of             genetree frequencies
  }
}
\details{
    this function will take a list of genetrees and call get_partial_cov_matrix to get the partial matrices from each genetree to        contstruct the final matrix
}
\value{
    returns an n x n matrix where n is the number of taxa in the genetree; uses genetree frequencies to scale covariances
}
\references{
    NEED TO ADD
}
\author{
    Lara Breithaupt and Mark Hibbins
}

}
\examples{
  final_genetree_matrix <- trees_to_vcv(test_input_files/seastar_genetrees_test_input.txt)
}

