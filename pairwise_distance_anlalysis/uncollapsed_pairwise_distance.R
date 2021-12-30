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
ps_tree <- read.tree("./ps_1524_uncollapsed_5_2018.nwk")
# extract the pairwise distance matrix from the ps_tree
ps_matrix <- cophenetic.phylo(ps_tree)

# Phage tree
phage_tree <- read.tree("./nj_MAFFT_cluster2_phages.nh")
# extract the pairise distance matrix from the phage_tree
phage_matrix <- cophenetic.phylo(phage_tree)



#################################################################################
############### cleaning up phage_matrix to look like ps_matrix #################
#################################################################################
colnames(phage_matrix) <- sub("late","", colnames(phage_matrix))
colnames(phage_matrix) <- sub("_pilon_contigs_renamed_phages_lysogenic","", colnames(phage_matrix))
colnames(phage_matrix) <- stri_replace_last_fixed(colnames(phage_matrix),'_','.')
colnames(phage_matrix) <- stri_replace_last_fixed(colnames(phage_matrix),'_',':')
colnames(phage_matrix) <- gsub(".*:","",colnames(phage_matrix))

# clean the colnames the same way
rownames(phage_matrix) <- sub("late","", rownames(phage_matrix))
rownames(phage_matrix) <- sub("_pilon_contigs_renamed_phages_lysogenic","", rownames(phage_matrix))
rownames(phage_matrix) <- stri_replace_last_fixed(rownames(phage_matrix),'_','.')
rownames(phage_matrix) <- stri_replace_last_fixed(rownames(phage_matrix),'_',':')
rownames(phage_matrix) <- gsub(".*:","",rownames(phage_matrix))



#################################################################################
###### delete ronames/ colnames in each df not present in the other df ##########
#################################################################################
# find all values in phage_tbl which are present in ps_matrix
phage_matrix2 <- phage_matrix[(colnames(phage_matrix)) %in% colnames(ps_matrix),(rownames(phage_matrix)) %in% rownames(ps_matrix)]

# find all values in ps_matrix which are present in phage_tbl2
ps_matrix2 <- ps_matrix[(colnames(ps_matrix)) %in% colnames(phage_matrix2),(rownames(ps_matrix)) %in% rownames(phage_matrix2)]



#################################################################################
################################# combine dfs ###################################
#################################################################################
# make a df where we have a column for phage and a column for pseeudomonas (instead of two matrices)
phage_melt <- melt(phage_matrix2)
colnames(phage_melt)[3] <- "phage"
ps_melt <- melt(ps_matrix2)
colnames(ps_melt)[3] <- "ps"
melt_both <- merge(phage_melt, ps_melt)



#################################################################################
###### plot Pseudomonas pairwise distance vs phage pairwise distance ############
#################################################################################
# make a column for labeling points in the next plot
melt_both$names <- paste(melt_both$Var1,melt_both$Var2)


# make a column to parse the genomes of interest (see where we are seeing clusters)
melt_both$group <- "not important"
melt_both$group[melt_both$Var1 %in% c("p26.F9","p25.A12","p26.D6","p25.D2")] <- "important"

cluster6 <- subset(melt_both, group == "important")


# plot pseudomonas pairwise distance vs phage pairwise distance
plot(melt_both$phage,melt_both$ps, col="lightblue", pch=19, cex=2)
text(melt_both$phage,melt_both$ps, labels=melt_both$names, cex=0.9, font=2)


library(ggplot2)
ggplot(melt_both, aes(x=phage, y=ps)) + geom_point() 




# subsetting portions of the phylogeny
# 1. Subsetting very divergent phage and very divergent ps (top right cluster of plot)