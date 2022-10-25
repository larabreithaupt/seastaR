remove(list=ls())

library("seastaR")
library("phytools")
# load trait matrix and newick string to tree
tomato_traits <- read.csv("/Users/larabreithaupt/Things/Hahn_Lab_Pruning_Algorithm/seastar/tomato_test/flower_morphometrics2.csv")
tomato_tree <- phytools::read.newick("/Users/larabreithaupt/Things/Hahn_Lab_Pruning_Algorithm//seastar/tomato_test/tomato_timetree.txt")

# convert from years to coales units
tomato_tree[["edge.length"]] <- tomato_tree[["edge.length"]] / 400000

# subset phylo object into high and low ILS knot groups:
# High ILS knot: S. galapagense 0436, S. cheesmaniae 3124, S. pimpinellifolium 1269
# Low ILS knot: S. pennellii 3778, S. pennellii 0716, S. pimpinellifolium 1589

low_ils_species <- c("S.galapagense", "S.cheesmaniae", "S.pim1269")
high_ils_species <- c("S.pen3778", "S.pen0716", "S.pim1589")
pruned_tree_low<-drop.tip(tomato_tree,tomato_tree$tip.label[-match(low_ils_species, tomato_tree$tip.label)])
pruned_tree_high<-drop.tip(tomato_tree,tomato_tree$tip.label[-match(high_ils_species, tomato_tree$tip.label)])

