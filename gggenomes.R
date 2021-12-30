library(tidyverse)
library(gggenomes)
library(ggplot2)
library(gggenes)
library(ggpubr)
library(ape)
library(ggtree)

# load in each subcluster of cluster 2 data
all <- read.table("./R_dat/protocol_and_files/tab_delimited_proteins/all_combined.csv", header = TRUE, sep="\t")

# make color palette for data with many classes
c25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)
pie(rep(1, 25), col = c25)

# plot all genes for all genomes using all df, but remove legend because it has too many values and won't plot
ggplot(all, aes(xmin = start, xmax = end, y = molecule, fill = gene)) +
  geom_gene_arrow() +
  facet_wrap(~ molecule, scales = "free", ncol = 1) +   theme(legend.position="none") 
# align genes across facets with make_alignment_dummies
dummies <- make_alignment_dummies(
  all,
  aes(xmin = start, xmax = end, y = molecule, id = gene),
  on = "Baseplate J-like protein"
)

ggplot(all, aes(xmin = start, xmax = end, y = molecule, fill = gene)) +
  geom_gene_arrow() +
  geom_blank(data = dummies) +
  facet_wrap(~ molecule, scales = "free", ncol = 1) + theme(legend.position = "none")

# subset to many dfs with similar sequences
#sub1
variabledata <- list("1_p4.G4_cl5","13_p6.B9_cl2","17_p20.D1_cl4","25_p25.G10_cl8","26_p20.C10_cl1")
x <- lapply(variabledata, function(x){subset(all, molecule == x)})
sub1 <- do.call(rbind, x)
#sub2
variabledata2 <- list("125_p25.C11_cl10","35_p20.H2_cl5","83_p13.E2_cl3")
x2 <- lapply(variabledata2, function(x2){subset(all, molecule == x2)})
sub2 <- do.call(rbind, x2)
#sub3
variabledata3 <- list("33_p2.H1_cl9","9_p13.D6_cl7")
x3 <- lapply(variabledata3, function(x3){subset(all, molecule == x3)})
sub3 <- do.call(rbind, x3)
#sub4
sub4 <- subset(all, molecule == "136_p25.A12_cl6")


# plot each subset
#sub1
p1 <- ggplot(sub1, aes(xmin = start, xmax = end, y = molecule, fill = gene)) +
  geom_gene_arrow() +
  facet_wrap(~ molecule, scales = "free", ncol = 1) + theme(legend.text = element_text(size=10)) +
  scale_fill_manual(values = c25)
p1
#sub2
p2 <- ggplot(sub2, aes(xmin = start, xmax = end, y = molecule, fill = gene)) +
  geom_gene_arrow() +
  facet_wrap(~ molecule, scales = "free", ncol = 1) + 
  scale_fill_manual(values = c25)
p2
#sub3
p3 <- ggplot(sub3, aes(xmin = start, xmax = end, y = molecule, fill = gene)) +
  geom_gene_arrow() +
  facet_wrap(~ molecule, scales = "free", ncol = 1) + 
  scale_fill_manual(values = c25)
p3
#sub4
p4 <- ggplot(sub4, aes(xmin = start, xmax = end, y = molecule, fill = gene)) +
  geom_gene_arrow() +
  facet_wrap(~ molecule, scales = "free", ncol = 1) + 
  scale_fill_manual(values = c25)
p4


# plot each subset aligned to a common gene
# sub1
table(sub1$gene)
dummies <- make_alignment_dummies(
  sub1,
  aes(xmin = start, xmax = end, y = molecule, id = gene),
  on = "Baseplate J-like protein"
)

ggplot(sub1, aes(xmin = start, xmax = end, y = molecule, fill = gene)) +
  geom_gene_arrow() +
  geom_blank(data = dummies) +
  facet_wrap(~ molecule, scales = "free", ncol = 1)

#sub2
table(sub2$gene)
dummies <- make_alignment_dummies(
  sub2,
  aes(xmin = start, xmax = end, y = molecule, id = gene),
  on = "Baseplate J-like protein"
)

ggplot(sub2, aes(xmin = start, xmax = end, y = molecule, fill = gene)) +
  geom_gene_arrow() +
  geom_blank(data = dummies) +
  facet_wrap(~ molecule, scales = "free", ncol = 1) + theme(legend.position = "none")

#sub3
table(sub3$gene)
dummies <- make_alignment_dummies(
  sub3,
  aes(xmin = start, xmax = end, y = molecule, id = gene),
  on = "Baseplate J-like protein"
)

ggplot(sub3, aes(xmin = start, xmax = end, y = molecule, fill = gene)) +
  geom_gene_arrow() +
  geom_blank(data = dummies) +
  facet_wrap(~ molecule, scales = "free", ncol = 1) + theme(legend.position = "none")

#sub4
table(sub4$gene)
dummies <- make_alignment_dummies(
  sub4,
  aes(xmin = start, xmax = end, y = molecule, id = gene),
  on = "Pyocin activator protein PrtN"
)

ggplot(sub4, aes(xmin = start, xmax = end, y = molecule, fill = gene)) +
  geom_gene_arrow() +
  geom_blank(data = dummies) +
  facet_wrap(~ molecule, scales = "free", ncol = 1) + theme(legend.position = "none")

# read in Ps tree 
tree <- read.tree("./R_dat/ps_1524_uncollapsed_5_2018.nwk")
# keep only the 10 subsetted genomes 
# make a df with only one column with all genomes we want to save
df <- read.csv("./R_dat/genomes.csv")
uneeded <- which(!(tree$tip.label) %in% df$genome_ID)
length(uneeded)
tree <- drop.tip(tree, uneeded)
plot(tree)

# see if we can plot the tree and gene plots together

