#install.packages('TreeDist')
library(TreeDist)
library(ape)
?TreeDist

# load in the pseudomonas OTU5 phylogeny 
ps_tree <- read.tree(file="./ps_1524_uncollapsed_5_2018.nwk")

# load in the core protein phylogeny you want to look at
PF04865.14_tree <- read.tree(file="./core_protein_trees/VOG01867.nh")

# clean tip labels to be the same
PF04865.14_tree$tip.label <- lapply(PF04865.14_tree$tip.label, gsub, pattern='_fa', replacement='')
PF04865.14_tree$tip.label <- gsub(".*p", "", PF04865.14_tree$tip.label)
PF04865.14_tree$tip.label <- paste0("p", PF04865.14_tree$tip.label)
PF04865.14_tree$tip.label[] <- lapply(PF04865.14_tree$tip.label, gsub, pattern='_', replacement='.')
PF04865.14_tree$tip.label <- as.character(PF04865.14_tree$tip.label)
# sanity check that they are the same now
#ps_tree$tip.label
#PF04865.14_tree$tip.label
# how many tips do our trees differ?
#length(PF04865.14_tree$tip.label) 
#length(ps_tree$tip.label) 

# drop tip labels from ps_tree that are not in core protein tree
uneeded <- which(!(ps_tree$tip.label) %in% PF04865.14_tree$tip.label) # which tip labels in core genome phylogeny are not in core protein tip labels
uneeded2 <- which(!(PF04865.14_tree$tip.label) %in% ps_tree$tip.label) # which tip labels in core protein tip labels are not in the core genome phylogeny
ps_tree2 <- drop.tip(ps_tree, uneeded)
coreprotein_tree <- drop.tip(PF04865.14_tree, uneeded2)
length(coreprotein_tree$tip.label)
length(ps_tree2$tip.label) # they are the same, yay!

# now use TreeDist
distance <- TreeDistance(coreprotein_tree, ps_tree2)
distance

# How similar are the two trees? -amount of phylogenetic information in common
SharedPhylogeneticInfo(coreprotein_tree, ps_tree2)
attr(SharedPhylogeneticInfo(coreprotein_tree, ps_tree2, reportMatching = TRUE), 'matching')
# this will be more useful with collapsed phylogeny
VisualizeMatching(SharedPhylogeneticInfo, coreprotein_tree, ps_tree2) # Which clades are matched?


# How different are the two trees?
DifferentPhylogeneticInfo(coreprotein_tree, ps_tree2) # Distance measure
DifferentPhylogeneticInfo(ps_tree2, coreprotein_tree) # The metric is symmetric
# both to see if it is symmetric

# Are they more similar than two trees of this shape would be by chance?
ExpectedVariation(coreprotein_tree, ps_tree2, sample=12)['DifferentPhylogeneticInfo', 'Estimate']


# Every split in tree1 conflicts with every split in tree3
# Pairs of conflicting splits contain clustering, but not phylogenetic, 
# information
SharedPhylogeneticInfo(coreprotein_tree, ps_tree2)
MutualClusteringInfo(coreprotein_tree, ps_tree2)

