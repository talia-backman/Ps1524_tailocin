# load libraries
library(ape)
library(ggplot2)
library(ggtree)
library(viridis)
library(tidyverse)

# read in data
# HTF haplotypes
hap_dat1 <- read.csv("./data/supplemental_data/HTF_aa_haplotypes.csv")
# TFA haplotypes
hap_dat2 <- read.csv("./data/supplemental_data/TFA_aa_haplotypes.csv")
# pseudomonas phylogeny
tree <- read.tree("./data/phylogeny_data/ps_1524_uncollapsed_5_2018.nwk")

# clean up hap_dat1
sort(unique(hap_dat1$HP12_haplotype),decreasing = TRUE)
# clean strain column
hap_dat1$Strain <- gsub(".*p", "", hap_dat1$Strain)
hap_dat1$Strain <- paste0("p", hap_dat1$Strain)

# merge data frames
hap_dat <- merge(hap_dat1, hap_dat2, by = "Strain")

# clean merged data frames for plotting
hap_dat$HP12_haplotype <- as.character(hap_dat$HP12_haplotype)
hap_dat$TFG_Haplotype <- as.character(hap_dat$TFG_Haplotype)

# drop tree tip labels that are not in hap_dat
uneeded <- which(!(tree$tip.label) %in% hap_dat$Strain)
tree <- drop.tip(tree, uneeded)

# remove duplicates based on Strain column
hap_dat <- hap_dat[!duplicated(hap_dat$Strain), ]
# make first column rownames in hap_dat
rownames(hap_dat) <- hap_dat[,1]
# remove remenant Strain column
hap_dat <- hap_dat[2:3]
# convert to matrix
hap_mat <- as.matrix(hap_dat)

# plot
p <- ggtree(tree, layout='circular', branch.length='none')
p3 <- gheatmap(p, hap_mat, colnames = FALSE, legend_title = "Cluster")
p_final <- p3 + scale_fill_viridis(discrete = TRUE) + 
  theme(legend.position = "left") + guides(fill=guide_legend(title= "TFG haplotype"))
p_final

# export to pdf
pdf("./figures/supplemental/TFA_HTF_MLphylogeny_linkage.pdf")
p_final
dev.off()