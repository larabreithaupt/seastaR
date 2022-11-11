remove(list=ls())
library(ape)
library(MASS)
library(phytools)
library(mvtnorm)
library(tidyverse)
library(geiger)
library(matlib)

### Functions ### 

sim_BM_fabio <- function(t1, t2, N, n_traits, sptree) {
 
 tau = t2 - t1
 scale = 2*N*0.0002 #Assuming constant mu
 
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

 var_covar <- scale*(unscaled_var_covar)
 
 traits <- mvrnorm(n = n_traits, mu = c(0, 0, 0), Sigma = var_covar)
 
 return(list(traits, rate_discordance, unscaled_var_covar))
 
}


ML_estimate_sigma2 <- function(trait, var_covar) {
  
  one <-c(1, 1, 1)
  zhat_root <- (t(one)%*%inv(var_covar)%*%one)%*%(t(one)%*%inv(var_covar)%*%trait)
  sigma2 <- ((t(trait - zhat_root*one))%*%inv(var_covar)%*%(trait - zhat_root*one))/3
  
  return(sigma2)
  
}


apply_ML_est <- function(dataset, right_var_covar) {

  right_ML_estimates <- apply(dataset, 1, ML_estimate_sigma2, right_var_covar)
  
  cov_matrix <- rep("Estimate", each = 1000)
  
  return(cbind(cov_matrix, right_ML_estimates))
}

### Analysis ### 

cond1_sptree <- ape::read.tree(
  "C:/Users/18126/OneDrive - University of Toronto/Projects/pruning_alg/sims_test/test_sptrees_sptreesims/cond1_sptree.txt")
cond2_sptree <- ape::read.tree(
  "C:/Users/18126/OneDrive - University of Toronto/Projects/pruning_alg/sims_test/test_sptrees_sptreesims/cond2_sptree.txt")
cond3_sptree <- ape::read.tree(
  "C:/Users/18126/OneDrive - University of Toronto/Projects/pruning_alg/sims_test/test_sptrees_sptreesims/cond3_sptree.txt")
cond4_sptree <- ape::read.tree(
  "C:/Users/18126/OneDrive - University of Toronto/Projects/pruning_alg/sims_test/test_sptrees_sptreesims/cond4_sptree.txt")
cond5_sptree <- ape::read.tree(
  "C:/Users/18126/OneDrive - University of Toronto/Projects/pruning_alg/sims_test/test_sptrees_sptreesims/cond5_sptree.txt")
cond6_sptree <- ape::read.tree(
  "C:/Users/18126/OneDrive - University of Toronto/Projects/pruning_alg/sims_test/test_sptrees_sptreesims/cond6_sptree.txt")
cond7_sptree <- ape::read.tree(
  "C:/Users/18126/OneDrive - University of Toronto/Projects/pruning_alg/sims_test/test_sptrees_sptreesims/cond7_sptree.txt")
cond8_sptree <- ape::read.tree(
  "C:/Users/18126/OneDrive - University of Toronto/Projects/pruning_alg/sims_test/test_sptrees_sptreesims/cond8_sptree.txt")
cond9_sptree <- ape::read.tree(
  "C:/Users/18126/OneDrive - University of Toronto/Projects/pruning_alg/sims_test/test_sptrees_sptreesims/cond9_sptree.txt")
cond10_sptree <- ape::read.tree(
  "C:/Users/18126/OneDrive - University of Toronto/Projects/pruning_alg/sims_test/test_sptrees_sptreesims/cond10_sptree.txt")
cond11_sptree <- ape::read.tree(
  "C:/Users/18126/OneDrive - University of Toronto/Projects/pruning_alg/sims_test/test_sptrees_sptreesims/cond11_sptree.txt")


condition1 <- sim_BM_fabio(4000, 50000, 2000, 1000, cond1_sptree)
condition2 <- sim_BM_fabio(4000, 50000, 6000, 1000, cond1_sptree)
condition3 <- sim_BM_fabio(4000, 50000, 8000, 1000, cond2_sptree)
condition4 <- sim_BM_fabio(4000, 50000, 12000, 1000, cond1_sptree)
condition5 <- sim_BM_fabio(4000, 50000, 16000, 1000, cond3_sptree)
condition6 <- sim_BM_fabio(4000, 50000, 24000, 1000, cond4_sptree)
condition7 <- sim_BM_fabio(4000, 50000, 32000, 1000, cond5_sptree)
condition8 <- sim_BM_fabio(4000, 50000, 40000, 1000, cond6_sptree)
condition9 <- sim_BM_fabio(4000, 50000, 48000, 1000, cond7_sptree)
condition10 <- sim_BM_fabio(4000, 50000, 80000, 1000, cond8_sptree)
condition11 <- sim_BM_fabio(4000, 50000, 150000, 1000, cond9_sptree)

condition1_results <- apply_ML_est(condition5[[1]], ape::vcv(cond5_sptree))
condition2_results <- apply_ML_est(condition5[[1]], condition2[[3]])
condition3_results <- apply_ML_est(condition5[[1]], condition3[[3]])
condition4_results <- apply_ML_est(condition5[[1]], condition4[[3]])
condition5_results <- apply_ML_est(condition5[[1]], condition5[[3]])
condition6_results <- apply_ML_est(condition5[[1]], condition6[[3]])
condition7_results <- apply_ML_est(condition5[[1]], condition7[[3]])
condition8_results <- apply_ML_est(condition5[[1]], condition8[[3]])
condition9_results <- apply_ML_est(condition5[[1]], condition9[[3]])
condition10_results <- apply_ML_est(condition5[[1]], condition10[[3]])
condition11_results <- apply_ML_est(condition5[[1]], condition11[[3]])

### Plot ### 

condition <- rep(c(0, 0.0144, 3.761e-02, 0.0981, 0.158, 0.255, 0.324, 0.375, 0.412, 
                   0.532, 0.592), 
                 each = 1000)

true_rate <- rep(2*16000*0.0002, each = 11000)

true <- rep("True", 11000)

all_estimates <- rbind(condition1_results, condition2_results, 
                       condition3_results, condition4_results,
                       condition5_results, condition6_results,
                       condition7_results, condition8_results,
                       condition9_results, condition10_results,
                       condition11_results)

truerates <- cbind(true, true_rate, condition)

results <- cbind(all_estimates, condition)
results <- rbind(results, truerates)

results_lineplot <- results %>%
  transform(estimates = as.numeric(right_ML_estimates),
            condition = as.numeric(condition),
            source = as.factor(cov_matrix)) %>%
  ggplot(aes(x=condition, y=estimates, group=source, color=source,
             linetype = source)) +
  stat_summary(fun.y=mean, geom="line", size=1.1) +
  scale_linetype_manual(values = c("solid", "dashed")) + 
  stat_summary(fun.data=mean_se, geom="pointrange") +
  geom_vline(xintercept = 0.158, size = 1) + 
  labs(x = "Rate of discordance used for estimate", y = "Rate estimate",
       color = "Source") + 
  theme_bw(base_size = 20) + 
  guides(linetype = "none")

results_lineplot
