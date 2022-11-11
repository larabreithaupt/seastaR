remove(list=ls())
library(ape)
library(MASS)
library(phytools)
library(mvtnorm)
library(tidyverse)
library(geiger)
library(matlib)

### Functions ### 

sim_BM_fabio <- function(t1, t2, N, n_traits) {
 
 tau = t2 - t1
 scale = 2*N*0.0002 #Assuming constant mu*sigma2
 
 rate_discordance = (2/3)*exp(-tau/(2*N))
 
 cov_AB = (1 - exp(-tau/(2*N)))*(1 + 
                (tau/(2*N) - (1 - (tau/(2*N))/(exp(tau/(2*N))-1)))) + 
                (1/3)*exp(-tau/(2*N))
 
 cov_nonsister = (1/3)*exp(-tau/(2*N))
 
 var = t1/(2*N) + (1 - exp(-tau/(2*N)))*(tau/(2*N) + 1) + 
            exp(-tau/(2*N))*(tau/(2*N) + 1 + 1/3)
 
 unscaled_var_covar <- matrix(c(var, cov_AB, cov_nonsister,
                                cov_AB, var, cov_nonsister,
                                cov_nonsister, cov_nonsister, var), 
                              nrow = 3, ncol = 3)
 
 var_covar <- scale*matrix(c(var, cov_AB, cov_nonsister,
                       cov_AB, var, cov_nonsister,
                       cov_nonsister, cov_nonsister, var), 
                     nrow = 3, ncol = 3)

 
 traits <- mvrnorm(n = n_traits, mu = c(0, 0, 0), Sigma = var_covar)
 
 return(list(traits, rate_discordance, unscaled_var_covar))
 
}


ML_estimate_sigma2 <- function(trait, var_covar) {
  
  one <-c(1, 1, 1)
  zhat_root <- (t(one)%*%inv(var_covar)%*%one)%*%(t(one)%*%inv(var_covar)%*%trait)
  sigma2 <- ((t(trait - zhat_root*one))%*%inv(var_covar)%*%(trait - zhat_root*one))/3
  
  return(sigma2)
  
}


apply_ML_est <- function(dataset, sptree, right_var_covar) {
  
  wrong_var_covar <- ape::vcv(sptree)

  wrong_ML_estimates <- apply(dataset, 1, ML_estimate_sigma2, wrong_var_covar)
  right_ML_estimates <- apply(dataset, 1, ML_estimate_sigma2, right_var_covar)
  
  
  cov_matrix <- rep(c("Species tree", "Gene tree"), each = 1000)
  estimates <- c(wrong_ML_estimates, right_ML_estimates)
  
  return(cbind(cov_matrix, estimates))
}

### Analysis ### 

condition1 <- sim_BM_fabio(4000, 50000, 2000, 1000)
print("Condition 1 simulations done")
condition2 <- sim_BM_fabio(4000, 50000, 4000, 1000)
print("Condition 2 simulations done")
condition3 <- sim_BM_fabio(4000, 50000, 6000, 1000)
print("condition 3 simulations done")
condition4 <- sim_BM_fabio(4000, 50000, 8000, 1000)
print("Condition 4 simulations done")
condition5 <- sim_BM_fabio(4000, 50000, 10000, 1000)
print("Condition 5 simulations done")
condition6 <- sim_BM_fabio(4000, 50000, 12000, 1000)
print("Condition 6 simulations done")
condition7 <- sim_BM_fabio(4000, 50000, 14000, 1000)
print("Condition 7 simulations done")

cond1_sptree <- ape::read.tree(
  "C:/Users/18126/OneDrive - University of Toronto/Projects/pruning_alg/sims_test/test_sptrees/cond1_sptree.txt")
cond2_sptree <- ape::read.tree(
  "C:/Users/18126/OneDrive - University of Toronto/Projects/pruning_alg/sims_test/test_sptrees/cond2_sptree.txt")
cond3_sptree <- ape::read.tree(
  "C:/Users/18126/OneDrive - University of Toronto/Projects/pruning_alg/sims_test/test_sptrees/cond3_sptree.txt")
cond4_sptree <- ape::read.tree(
  "C:/Users/18126/OneDrive - University of Toronto/Projects/pruning_alg/sims_test/test_sptrees/cond4_sptree.txt")
cond5_sptree <- ape::read.tree(
  "C:/Users/18126/OneDrive - University of Toronto/Projects/pruning_alg/sims_test/test_sptrees/cond5_sptree.txt")
cond6_sptree <- ape::read.tree(
  "C:/Users/18126/OneDrive - University of Toronto/Projects/pruning_alg/sims_test/test_sptrees/cond6_sptree.txt")
cond7_sptree <- ape::read.tree(
  "C:/Users/18126/OneDrive - University of Toronto/Projects/pruning_alg/sims_test/test_sptrees/cond7_sptree.txt")


condition1_results <- apply_ML_est(condition1[[1]], cond1_sptree, condition1[[3]])
print("Condition 1 analyses done")
condition2_results <- apply_ML_est(condition2[[1]], cond2_sptree, condition2[[3]])
print("Condition 2 analyses done")
condition3_results <- apply_ML_est(condition3[[1]], cond3_sptree, condition3[[3]])
print("Condition 3 analyses done")
condition4_results <- apply_ML_est(condition4[[1]], cond4_sptree, condition4[[3]])
print("Condition 4 analyses done")
condition5_results <- apply_ML_est(condition5[[1]], cond5_sptree, condition5[[3]])
print("Condition 5 analyses done")
condition6_results <- apply_ML_est(condition6[[1]], cond6_sptree, condition6[[3]])
print("Condition 6 analyses done")
condition7_results <- apply_ML_est(condition7[[1]], cond7_sptree, condition7[[3]])
print("COndition 7 analysis done")

### Plot ### 

condition <- rep(c(condition1[[2]], condition2[[2]],
                   condition3[[2]], condition4[[2]],
                   condition5[[2]], condition6[[2]],
                   condition7[[2]]), 
                   each = 2000)

true_rate <- rep(c(2*2000*0.0002, 2*4000*0.0002, 2*6000*0.0002,
                   2*8000*0.0002, 2*10000*0.0002, 2*12000*0.0002,
                   2*14000*0.0002), each = 2000)

true <- rep("True", 14000)

all_estimates <- rbind(condition1_results, condition2_results, 
                       condition3_results, condition4_results,
                       condition5_results, condition6_results,
                       condition7_results)

truerates <- cbind(true, true_rate, condition)

results <- cbind(all_estimates, condition)
results <- rbind(results, truerates)

results_lineplot <- results %>%
  transform(estimates = as.numeric(estimates),
            condition = as.numeric(condition),
            source = as.factor(cov_matrix)) %>%
  ggplot(aes(x=condition, y=estimates, group=source, color=source,
             linetype = source)) +
  stat_summary(fun.y=mean, geom="line", size=1.1) +
  scale_linetype_manual(values = c("solid", "solid", "dashed")) + 
  stat_summary(fun.data=mean_se, geom="pointrange") +
  labs(x = "Rate of discordance", y = "Rate estimate",
       color = "Source") + 
  theme_bw(base_size = 20) + 
  guides(linetype = "none")

results_lineplot
