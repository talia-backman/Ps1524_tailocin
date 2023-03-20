# tutorial found here: https://evolutionarygenetics.github.io/Chapter7.html
# load libraries
library(tidyverse)
library(ape)
library(pegas)
library(spider)


# read in data
# dna seq
mydna <- read.dna("./data/supplemental_data/all_concat_1524.fasta", format = "fasta")

# create sliding windows
mydna_windows <- slidingWindow(mydna, width = 1000, interval = 100)

# how to extract one window and calculate tajima's d on it
#tajima.test(mydna_windows[[1]])

# loop through windows and calculate td
# then append to td_df
td_df = data.frame()
for(i in mydna_windows){
  stats <- data.frame(tajima.test(i))
  td_df <- rbind(td_df,stats)
}
# add window column
td_df$window <- (1:nrow(td_df))*100
# write to csv 
write.csv(td_df, "./data/supplemental_data/800bp_windows_allsites_td_dat.csv")

# now do the same thing but for pi
# loop through windows and calculate pi
# then append to pi_df
pi_df = data.frame()
for(i in mydna_windows){
  stats <- data.frame(nuc.div(i))
  pi_df <- rbind(pi_df,stats)
}
# add window column
pi_df$window <- (1:nrow(pi_df))*100
# write to csv 
write.csv(pi_df, "./data/supplemental_data/800bp_windows_allsites_pi_dat.csv")