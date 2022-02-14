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
    Triplet pulled from species tree input
}
  \item{sptree}{
      Species tree phylo object produced from file parser (parse_input_file); species tree is in coalescent units
  }
}
\details{
    this function will be called on each triplet to get the internal branch from each triplet in the species tree
}
\value{
    returns a list containing the length of the internal branch, the final height (the speceies tree height - the triplet height) and a string showing the sister taxa from the triplet
}
\references{
    NEED TO ADD
}
\author{
    Lara Breithaupt and Mark Hibbins
}

}
\examples{
  internal_branch <- get_internal_branch(triplets[, i], sptree)
}

