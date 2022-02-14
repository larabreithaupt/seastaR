\name{parse_input_file}
\alias{parse_input_file}

\title{
parse_input_file
}
\description{
    Function to parse the input file.
}
\usage{
    parse_input_file(file_path, genetrees=TRUE)
}

\arguments{
  \item{file_path}{
      Path to the input file; ensure that last line of file has line break
}
  \item{genetrees}{
      Boolean flag to specify genetree vs species tree input; defaults to true
  }
}
\details{
    Species tree input is a single line Newick string and genetree input is a Nexus file input with the genetrees and the genetree       frequencies under seperate headers
}
\value{
    If genetrees is false, it returns a single phylo object and if true, it returns a list of phylo objects where the last element of     the list is a vector of the genetree frequencies
}
\references{
        Paradis E, Schliep K (2019). “ape 5.0: an environment for modern phylogenetics and evolutionary analyses in R.” Bioinformatics, 35, 526-528.
}
\author{
    Lara Breithaupt and Mark Hibbins
}

}
\examples{
sptree_phylo_object <- parse_input_file(test_input_files/seastar_sptree_test_input.txt, genetrees = FALSE)

genetree_phylo_object_list <- parse_input_file(test_input_files/seastar_genetrees_test_input)

}

