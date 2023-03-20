# load libraries
library(ggtree)
library(ggplot2)
library(viridis)
library(ape)

# read in data and clean it up for plotting
df <- read.csv("./data/mash_data/clustering_dat_clean_ordered.csv", header = F)
all_df_matrix <- as.matrix(df)
rownames(all_df_matrix) <- all_df_matrix[,1]
all_df_matrix <- all_df_matrix[ ,c(2:10)]

# read in tree data (pseudomonas core genome maximum likelihood phylogeny) and clean it up 
tree <- read.tree(file="./data/phylogeny_data/ps_1524_uncollapsed_5_2018.nwk")
# cleaning the tree: remove rows in tree that aren't in df1 
uneeded <- which(!(tree$tip.label) %in% df$V1)
length(uneeded) # there are 11 genomes with no viral elements in them
tree <- drop.tip(tree, uneeded)

# plot clusters to phylogeny
p <- ggtree(tree, layout='circular', branch.length='none')
p3 <- gheatmap(p, all_df_matrix, colnames = FALSE, legend_title = "Cluster") + scale_fill_viridis(discrete = TRUE)
p4 <- p3 + geom_strip('p27.D6','p22.D1', fontsize = 7, barsize=1, offset=75) +
  geom_strip('p2.B6','p6.B5', fontsize = 7,barsize = 1, offset = 75) + 
  theme(legend.position = "none")

# export figure as a pdf
pdf("./figures/01_phage_annotation/viral_clusters_to_phylogeny.pdf")
p4
dev.off()