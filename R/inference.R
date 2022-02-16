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
