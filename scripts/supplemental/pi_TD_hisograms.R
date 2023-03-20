# load libraries
library(ggplot2)

# read in data
dat <- read.csv("./data/supplemental_data/800bp_windows_allsites_pi_dat.csv")
dat2 <- read.csv("./data/supplemental_data/800bp_windows_allsites_td_dat.csv")

# plot and export as pdf
pdf("./figures/supplemental/pi_histogram.pdf")
ggplot(dat, aes(x=segsites)) + geom_histogram(fill = "#440154FF") + 
  xlab("Nucleotide diversity") +
  ylab("Frequency")
dev.off()

# plot and export as pdf
pdf("./figures/supplemental/TD_histogram.pdf")
ggplot(dat2, aes(x=D)) + geom_histogram(fill = "#440154FF") + 
  xlab("Tajima's D") +
  ylab("Frequency")
dev.off()