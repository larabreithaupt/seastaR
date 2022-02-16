\name{sigma2_inference}
\alias{sigma2_inference}

\title{
sigma2_inference
}
\description{
      Takes a trait and a covariance matrix and returns a maximum likelihood estimate of sigma squared
}
\usage{
    sigma2_inference(var_covar, trait)
}

\arguments{
  \item{var_covar}{
    variance, covariance matrix that is #taxa x #taxa in size
  }
  \item{trait}{
    a numeric vector of size #taxa that include the trait values of each taxon; the order of the values in the vector must match the order of the rows in the variance, covariance matrix
  }

}
\details{
    This function will return an estimated sigma squared value based on the un-biased maximum likelihood estimator as described in O'Meara et al. (2006).
}
\value{
    Returns the maximum likelihood estimate of sigma squared
}
\references{
    O'Meara BC, Ane C, Sanderson MJ, Wainwright PC. 2006. Testing for different rates of continuous trait evolution using likelihood. Evolution 60(5): 922-933
}
\author{
    Lara Breithaupt and Mark Hibbins
}
\examples{
  sptree <- parse_input_file(test_input_files/seastar_sptree_test_input.txt, genetrees = FALSE)
  var_covar <- get_full_matrix(sptree)

  traits <- c(1, 2, 3)

  estimated_sigma2 <- sigma2_inference(var_covar, traits)
}

