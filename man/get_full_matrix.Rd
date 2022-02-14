\name{get_full_matrix}
\alias{get_full_matrix}

\title{
get_full_matrix
}
\description{
    This function returns the final species tree matrix
}
\usage{
    get_full_matrix(sptree)
}

\arguments{
  \item{sptree}{
     Species tree phylo object produced from file parser (parse_input_file); species tree is in coalescent units
}
\value{
  returns a final nxn matrix containing variance and covariance values
}
\references{
  NEED TO ADD
}
\author{
    Lara Breithaupt and Mark Hibbins
}

}
\examples{
  full_matrix <- get_full_matrix(sptree)
}
