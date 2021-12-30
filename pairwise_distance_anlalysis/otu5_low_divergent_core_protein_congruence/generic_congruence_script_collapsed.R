library(ape)
library(dplyr)
library(treeio)
library(castor)

# load in the pseudomonas OTU5 phylogeny 
ps_tree <- read.tree(file="./collapsed_single_name.nwk")

# load in the core protein phylogeny you want to look at
PF04865.14_tree <- read.tree(file="./core_protein_trees/VOG01059.nh")

length(ps_tree$tip.label)
length(PF04865.14_tree$tip.label)

# clean tip labels to be the same then drop tip labels from ps_tree that are not in core protein tree
PF04865.14_tree$tip.label <- lapply(PF04865.14_tree$tip.label, gsub, pattern='_fa', replacement='')
PF04865.14_tree$tip.label <- gsub(".*p", "", PF04865.14_tree$tip.label)
PF04865.14_tree$tip.label <- paste0("p", PF04865.14_tree$tip.label)
PF04865.14_tree$tip.label[] <- lapply(PF04865.14_tree$tip.label, gsub, pattern='_', replacement='.')
PF04865.14_tree$tip.label <- as.character(PF04865.14_tree$tip.label)
ps_tree$tip.label
PF04865.14_tree$tip.label

length(PF04865.14_tree$tip.label) 
length(ps_tree$tip.label) 

uneeded <- which(!(ps_tree$tip.label) %in% PF04865.14_tree$tip.label) # which tip labels in core genome phylogeny are not in core protein tip labels
uneeded2 <- which(!(PF04865.14_tree$tip.label) %in% ps_tree$tip.label) # which tip labels in core protein tip labels are not in the core genome phylogeny
ps_tree2 <- drop.tip(ps_tree, uneeded)
coreprotein_tree <- drop.tip(PF04865.14_tree, uneeded2)
length(coreprotein_tree$tip.label)
length(ps_tree2$tip.label) 

# make the 'association' dataframe 
df <- coreprotein_tree$tip.label
df <- as.data.frame(df)
# rename first column to tip_labels
df <- df %>% rename(
  tip_labels = df)
# create a new column with duplicated data
df$associations <- df$tip_labels

# plot two phylogenetic trees face to face with links between
cophyloplot(ps_tree2, coreprotein_tree, assoc = df, space=250, gap=5, col='blue', show.tip.label = FALSE)
