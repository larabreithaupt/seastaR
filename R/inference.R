#import(MASS)

sigma2_inference <- function(var_covar, trait){

  num_taxa <- length(var_covar[, 1])
  if (num_taxa == length(trait)){

    one <- c(rep(1, num_taxa))

    zhat_root <- (t(one)%*%solve(var_covar)%*%one)%*%(t(one)%*%solve(var_covar)%*%trait)
    sigma2 <- ((t(trait - zhat_root*one))%*%solve(var_covar)%*%(trait - zhat_root*one))/(num_taxa-1)

    return(sigma2)
  }
  else{
    stop("ERROR, trait vector is not the length of the number of taxa")
  }

}

simulate_traits <- function(n_traits, var_covar, sigma2){

  num_taxa <- length(var_covar[, 1])
  traits <- mvrnorm(n = n_traits, mu = c(rep(0, num_taxa)), Sigma = sigma2*var_covar)
  return(traits)
}


sigma2_likelihood <- function(sigma2, trait, var_covar) {

  one <-c(1, 1, 1)
  zhat_root <- (t(one)%*%solve(var_covar)%*%one)%*%(t(one)%*%solve(var_covar)%*%trait)
  top <- exp((-1/2)*((t(trait - zhat_root*one))%*%solve(sigma2*var_covar)%*%(trait - zhat_root*one)))
  bottom = sqrt(((2*pi)^3)*det(sigma2*var_covar))
  logL = log(top/bottom)

  return(logL)
}

trait_likelihood_surface <- function(trait, var_covar, lower_sigma2, higher_sigma2){

  trait_sigma2_vals <- seq(lower_sigma2, higher_sigma2, length.out = 1000)
  likelihoods <- lapply(trait_sigma2_vals, sigma2_likelihood, trait, var_covar)

  return(unlist(likelihoods))

}

