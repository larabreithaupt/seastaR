remove(list=ls())

library("seastaR")
library("phytools")
library("tidyverse")
library("ape")

# load trait matrix and newick string to tree
tomato_traits <- read.csv("/Users/larabreithaupt/Things/Hahn_Lab_Pruning_Algorithm/seastar/tomato_test/flower_morphometrics2.csv")
tomato_tree <- phytools::read.newick("/Users/larabreithaupt/Things/Hahn_Lab_Pruning_Algorithm//seastar/tomato_test/tomato_timetree.txt")

# convert from years to coales units
tomato_tree[["edge.length"]] <- tomato_tree[["edge.length"]] / 400000

# subset phylo object into high and low ILS knot groups:
# High ILS knot: S. galapagense 0436, S. cheesmaniae 3124, S. pimpinellifolium 1269
# Low ILS knot: S. pennellii 3778, S. pennellii 0716, S. pimpinellifolium 1589

high_ils_species <- c("S.galapagense", "S.cheesmaniae", "S.pim1269")
low_ils_species <- c("S.pen3778", "S.pen0716", "S.pim1589")
pruned_tree_low<-drop.tip(tomato_tree,tomato_tree$tip.label[-match(low_ils_species, tomato_tree$tip.label)])
pruned_tree_high<-drop.tip(tomato_tree,tomato_tree$tip.label[-match(high_ils_species, tomato_tree$tip.label)])

# get species tree covar matrix
low_ils_C <- ape::vcv(pruned_tree_low)
high_ils_C <- ape::vcv(pruned_tree_high)

#format traits and get full var covar matrices
high_ils_Cstar <- get_full_matrix(pruned_tree_high)
low_ils_Cstar <- get_full_matrix(pruned_tree_low)

low_ils_accession <- c("LA1589", "LA3778", "LA0716")
high_ils_accession <- c("LA0436", "LA3124", "LA1269")

low_ils_filtered_traits <- dplyr::filter(tomato_traits, AccessionID %in% low_ils_accession)
high_ils_filtered_traits <- dplyr::filter(tomato_traits, AccessionID %in% high_ils_accession)

# calculate mean within each accession
low_ils_mean_acc <- low_ils_filtered_traits %>%
  dplyr::group_by(AccessionID) %>%
  summarize_at(vars(SE_PlantMean:SE_bin), list(name = mean))
high_ils_mean_acc <- high_ils_filtered_traits %>%
  dplyr::group_by(AccessionID) %>%
  summarize_at(vars(SE_PlantMean:SE_bin), list(name = mean))

# get specific traits
low_ils_traits <- low_ils_mean_acc[ ,c(1, 3:5)]
high_ils_traits <- high_ils_mean_acc[ ,c(1, 3:5)]

# reorder rows of trait df to match covar matrix
low_ils_traits <- low_ils_traits %>% dplyr::slice(match(low_ils_accession, AccessionID))
high_ils_traits <- high_ils_traits %>% dplyr::slice(match(high_ils_accession, AccessionID))

# estimate sigma^2 value

low_ils_AL_C_s2 <- sigma2_inference(low_ils_C, unlist(low_ils_traits[, 3]))





