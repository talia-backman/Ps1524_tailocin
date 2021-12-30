library(ape)
library(dplyr)
library(treeio)
library(castor)

# load in the pseudomonas OTU5 phylogeny 
ps_tree <- read.tree(file="./ps_1524_uncollapsed_5_2018.nwk")

# load in the core protein phylogeny you want to look at
PF04865.14_tree <- read.tree(file="./core_protein_trees/VOG10943.nh")

# clean tip labels to be the same then drop tip labels from ps_tree that are not in core protein tree
PF04865.14_tree$tip.label <- lapply(PF04865.14_tree$tip.label, gsub, pattern='_fa', replacement='')
PF04865.14_tree$tip.label <- gsub(".*p", "", PF04865.14_tree$tip.label)
PF04865.14_tree$tip.label <- paste0("p", PF04865.14_tree$tip.label)
PF04865.14_tree$tip.label[] <- lapply(PF04865.14_tree$tip.label, gsub, pattern='_', replacement='.')
PF04865.14_tree$tip.label <- as.character(PF04865.14_tree$tip.label)
#ps_tree$tip.label
#PF04865.14_tree$tip.label

length(PF04865.14_tree$tip.label) # 1398
length(ps_tree$tip.label) # 1524

uneeded <- which(!(ps_tree$tip.label) %in% PF04865.14_tree$tip.label)
length(uneeded)
ps_tree2 <- drop.tip(ps_tree, uneeded)
length(ps_tree2$tip.label) # 1398

# make the 'association' dataframe 
df <- PF04865.14_tree$tip.label
df <- as.data.frame(df)
# rename first column to tip_labels
df <- df %>% rename(
  tip_labels = df)
# create a new column with duplicated data
df$associations <- df$tip_labels

# plot two phylogenetic trees face to face with links between
cophyloplot(ps_tree2, PF04865.14_tree, assoc = df, space=4000, gap=50, col='blue', show.tip.label = FALSE)


# random subset to links in order to observe patterns better
set.seed(123)
test <- sample_n(df, size = 100)
cophyloplot(ps_tree2, PF04865.14_tree, assoc = test, length.line = 50, space=4000, gap=50, col='blue', show.tip.label = FALSE)

