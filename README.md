# seastaR

Source code for the R package `seastaR`, written by Lara Breithaupt and Mark Hibbins (see Hibbins et al. 2022). This package builds a gene tree variance-covariance matrix that can be used for downstream phylogenetic comparative analyses. Also included is a trait simulator and a function for maximum-likelihood estimation of the evolutionary rate.

## Installation 

Install from GitHub with `devtools::install_github("larabreithaupt/seastaR")`.

## Usage example

See below for a basic usage example highlighting the main functions. First we load in a species phylogeny in coalescent units:

    > library(seastaR)
    > library(ape)
    > sptree_example <- seastaR::parse_input_file("tests/seastaR_sptree_example.txt", genetrees = FALSE)
    
Next we construct both a standard species tree covariance matrix using ape, and a gene tree covariance matrix using our `get_full_matrix()` approach, which breaks the tree down into triplets and uses theory to build the matrix:

    > C_matrix <- ape::vcv(sptree_example) #Standard species tree matrix 
    > Cstar_matrix <- seastaR::get_full_matrix(sptree_example) #Gene tree matrix 

We can simulate a test trait from the gene tree covariance matrix:

    > test_trait <- seastaR::simulate_traits(1, Cstar_matrix, 1)

Finally, we use our implementation of the maximum-likelihood estimator from O'Meara et al. 2006 to estimate the evolutionary rate:

    > sptree_rate <- seastaR::sigma2_inference(C_matrix, test_trait) #Standard estimate
    > sptree_rate
             [,1]
    [1,] 1.835745
    > genetree_rate <- seastaR::sigma2_inference(Cstar_matrix, test_trait) #Gene tree estimate
    > genetree_rate
             [,1]
    [1,] 1.188734

We can do a similar analyses using a specified set of input gene trees instead. First we read the input file:

    > genetree_example <- seastaR::parse_input_file("tests/seastar_genetrees_example.txt", genetrees = TRUE)

Then we use `trees_to_vcv()` to construct the gene tree covariance matrix by weighting all the branches of the input gene trees: 

    > Cstar_matrix <- seastaR::trees_to_vcv(genetree_example)
    > Cstar_matrix 
             sp2      sp3      sp4
    sp2 2.382268 0.183000 0.183000
    sp3 0.183000 2.382268 0.782830
    sp4 0.183000 0.782830 2.382268

Simulate a test trait:

    > test_trait <- seastaR::simulate_traits(1, Cstar_matrix, 1)

And estimate the evolutionary rate: 

    > #Rate estimate 
    > genetree_rate <- seastaR::sigma2_inference(Cstar_matrix, test_trait) #Gene tree estimate
    > genetree_rate
              [,1]
    [1,] 0.5903072

## Analyses

The `analyses/` folder contains several scripts used for analyses in Hibbins et al. 2022. `accuracy_ML_est.R` tests the accuracy of the gene tree covariance matrix method using simulations and generated the results in Figure 4A. `gt_error_ML_est.R` tested the effects of gene tree estimation error on the covariance matrix method using simulations, and generated the results in Figure 5A. `tomato_seastar_analysis.R` uses seastaR to estimate evolutionary rates for empirical tomato floral traits, and generated the results in Figure 6A. 
