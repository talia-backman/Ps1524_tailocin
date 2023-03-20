# load libraries
library(ggplot2)
library(gggenes) 
library(ggpubr)
library(ape)
library(ggtree)

# make color palette 
c25 <- c("#440154FF", "#482677FF","#404788FF","#33638DFF","#2D708EFF","#238A8DFF","#29AF7FFF","#55C667FF","#73D055FF","#B8DE29FF","#DCE319FF","#FDE725FF","#f0f397","#e2e82e","#d8b427","#36aa6a","#1a7f75","grey","black","#378aae","#8f81c9","#680280","#525aad","#3a2e6c")

# read in data
# tab delimited gene data
all <- read.table("./data/gene_data/ps41_gene_data.csv", header = TRUE, sep=",")
# load in df that has a list of the genome IDs used (to subset tree with)
df <- read.csv("./data/gene_data/41genomes.csv")
# Pseudomonas phylogeny 
tree <- read.tree("./data/phylogeny_data/ps_1524_uncollapsed_5_2018.nwk")
# drop the tree$tip.label which are not in df$genome_ID
uneeded <- which(!(tree$tip.label) %in% df$genome_ID)
length(uneeded)
tree <- drop.tip(tree, uneeded)

# plot the tree
p_tree <- ggtree(tree) + geom_tiplab() + geom_treescale()
# with clade labels
p_tree2 <- p_tree + 
  geom_strip('p25.A12','p8.H7', fontsize = 5, barsize=1,
             label='ATUE5',offset=0.5, offset.text=0.1,angle = 90) +
  geom_strip('p6.G2','p6.B5', fontsize = 5, barsize=1,
             label='Other OTUs',offset=0.5,offset.text=0.1, angle = 90,hjust =.5)

# gene plot
p1 <- ggplot(all, aes(xmin = start, xmax = end, y = molecule, fill = gene)) +
  geom_gene_arrow() + scale_fill_manual(values = c25) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.line.y =element_blank(), axis.line.x = element_blank(), axis.ticks = element_blank(),
        axis.text=element_blank(),
        axis.title.y=element_blank(),
        panel.border=element_blank(),
        plot.background=element_blank()) 

# combine tree and gene plots into one figure
# export as pdf
pdf("./figures/02_gene_content_VC2/gene_plot.pdf")
ggarrange(p_tree2, p1, ncol = 2, nrow = 1, widths = c(4,7))
dev.off()
