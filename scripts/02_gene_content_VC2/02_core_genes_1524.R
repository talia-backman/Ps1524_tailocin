# load libraries
library(ape)
library(ggplot2)
library(dplyr)

# read in data
# pseudomonas phylogeny
tree <- read.tree("./data/phylogeny_data/ps_1524_uncollapsed_5_2018.nwk")
# gene presence absence table (1 =  ortholog present, 0 = ortholog absent)
dat <- read.csv("./data/gene_data/gene_presence_absence1524.csv")

# clean data
# make ortholog ID the rownames
rownames(dat) <- dat$ortholog_id
# remove ortholog_id column from dat
dat <- dat[-c(1)]
# also remove annot column from dat
dat <- dat[-c(1525)]

# calculate absolute totals and percent totals
# make a column at the end of dat which adds how many times each ortholog is present
dat_totals <- cbind(dat, Total = rowSums(dat))
# make a column at the end of dat which adds percent each ortholog is present
dat_percents <- cbind(dat, Percentage = ((rowSums(dat))/1524)*100)

# plot a histogram of frequency of the orthologs in these genomes (with absolute totals)
p1 <- ggplot(dat_totals, aes(x=Total)) + geom_histogram(color = "black",fill = "#404788FF") + 
  scale_y_continuous(trans='log10') + theme_classic() +
  labs(x = "Number of Genomes", y = "Number of Genes (log10)") 
# plot a histogram of frequency of the orthologs in these genomes (with percent totals)
p2 <- ggplot(dat_percents, aes(x=Percentage)) + geom_histogram(color = "black",fill = "#404788FF") + 
  scale_y_continuous(trans='log10') + theme_classic() +
  labs(x = "Percent of Genomes", y = "Number of Genes (log10)")

# export as pdfs
pdf("./figures/02_gene_content_VC2/1524_pan_genome_absolute_totals.pdf")
p1
dev.off()
pdf("./figures/02_gene_content_VC2/1524_pan_genome_percent_totals.pdf")
p2
dev.off()
