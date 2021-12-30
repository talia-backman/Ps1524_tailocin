# load libraries
library(ape)
library(castor)
library(ggplot2)
library(ggtree)
library(dplyr)
library(flexclust)
library(quanteda)
library(spaa)
library(reshape2)
library(stringi)


#################################################################################
######## getting pairwise distances of cluster 2 phage and Pseudomonas ##########
#################################################################################

# collapsed Pseudomonas tree
ps_tree <- read.tree("./collapsed_single_name.nwk")
ps_tree
plot(ps_tree)
length(ps_tree$tip.label) # 165
subset <- ps_tree$tip.label
subset <- as.data.frame(subset)
write.csv(subset, "./collapsed_tree_tip_labels.csv")

otu5 <- read.csv("./low_divergent_cluster_otu5.csv")

# subset tree to only otu5 strains
uneeded <- which(!(ps_tree$tip.label) %in% otu5$genome_ID) # this is actually needed
tree <- drop.tip(ps_tree, uneeded)
plot(tree)
length(tree$tip.label)

otu5_subset <- tree$tip.label
otu5_subset <- as.data.frame(otu5_subset)
write.csv(otu5_subset, "./collapsed_tree_otu5.csv")

write.tree(tree, "./collapsed_otu5.nwk")


