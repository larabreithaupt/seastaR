\name{get_triplet_branches}
\alias{get_triplet_branches}

\title{
get_triplet_branches
}
\description{
    This function returns a data frame containing the internal branch lengths for all triplets contained in the species tree
}
\usage{
    get_triplet_branches(sptree)
}

\arguments{
  \item{sptree}{
     Species tree phylo object produced from file parser (parse_input_file); species tree is in coalescent units
  }
}
\value{
    returns a data frame with the triplet name, its tau value, and its internal branches' height
}
\references{
    NEED TO ADD
}
\author{
    Lara Breithaupt and Mark Hibbins
}

\examples{
    branches_df <- get_triplet_branches(sptree)
}

