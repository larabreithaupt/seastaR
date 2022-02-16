\name{simulate_traits}
\alias{simulate_traits}

\title{
simulate_traits
}
\description{
    Simulates trait values from an input matrix and sigma squared
}
\usage{
    simulate_traits(n_traits, var_covar, sigma2)
}

\arguments{
  \item{n_traits}{
    integer value for the number of traits to simulate
  }
  \item{var_covar}{
    variance, covariance matrix that is #taxa x #taxa in size
  }
  \item{sigma2}{
    sigma squared value
  }
}
\details{
    This function simulates traits by drawing from a multivariate normal distribution
}
\value{
    Returns a matrix of size #traits x #taxa with trait values as entries
}
\author{
    Lara Breithaupt and Mark Hibbins
}

\examples{

simulates <- simulate_traits(4, matrix(c(1,1,1,
                                       1,1,1,
                                       1,1,1), nrow = 3, ncol = 3), .5)
}

