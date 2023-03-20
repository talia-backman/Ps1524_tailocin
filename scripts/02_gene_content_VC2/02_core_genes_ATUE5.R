# load libraries
library(ape)
library(ggplot2)
library(dplyr)

# read in data
# pseudomonas phylogeny
tree <- read.tree("./data/phylogeny_data/ps_1524_uncollapsed_5_2018.nwk")
# gene presence absence table (1 =  ortholog present, 0 = ortholog absent)
dat <- read.csv("./data/gene_data/gene_presence_absence1524.csv")
# ATUE IDs data
ATUEs <- read.csv("./data/gene_data/ATUE_IDs.csv", sep = "\t")

# subset to only ATUE5 (or ATUE1 in this dataframe)
ATUE1 <- subset(ATUEs, Otu == "1")
# remove rows in dat that aren't in ATUE1 
uneeded <- which(!(colnames(dat)) %in% ATUE1$Strain)
dat_ATUE1 <- dat[-c(uneeded)]

# calculate absolute totals and percent totals
# make a column at the end of dat which adds how many times each ortholog is present
dat_totals <- cbind(dat_ATUE1, Total = rowSums(dat_ATUE1))
# make a column at the end of dat which adds percent each ortholog is present
dat_percents <- cbind(dat_ATUE1, Percentage = ((rowSums(dat_ATUE1))/1524)*100)

# plot a histogram of frequency of the orthologs in these genomes (with absolute totals)
p1 <- ggplot(dat_totals, aes(x=Total)) + geom_histogram(color = "black",fill = "#404788FF") + 
  scale_y_continuous(trans='log10') + theme_classic() +
  labs(x = "Number of Genomes", y = "Number of Genes (log10)") 
# plot a histogram of frequency of the orthologs in these genomes (with percent totals)
p2 <- ggplot(dat_percents, aes(x=Percentage)) + geom_histogram(color = "black",fill = "#404788FF") + 
  scale_y_continuous(trans='log10') + theme_classic() +
  labs(x = "Percent of Genomes", y = "Number of Genes (log10)")

# export as pdfs
pdf("./figures/02_gene_content_VC2/ATUE5_pan_genome_absolute_totals.pdf")
p1
dev.off()
pdf("./figures/02_gene_content_VC2/ATUE5_pan_genome_percent_totals.pdf")
p2
dev.off()