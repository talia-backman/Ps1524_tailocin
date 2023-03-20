# tutorial found here: https://evolutionarygenetics.github.io/Chapter7.html
library(tidyverse)
library(ape)
library(pegas)
library(PopGenome)
library(msa)

# get a list of files
getwd()
setwd("/home/talia-backman/Desktop/Ps1524_tailocin/data/panX_data")
input_files <- list.files("/home/talia-backman/Desktop/Ps1524_tailocin/data/panX_data", pattern = ".pal2nal")

# loop for reading input and writing an output with population genetics summary statistics
df_total = data.frame()
for(i in 1:length(input_files)){
  file <- (input_files[i])
  filename <- gsub(".pal2nal", "", file)
  mydna <- read.dna(file, format = "fasta")
  pi <- nuc.div(mydna)
  ss <- length(seg.sites(mydna))
  S <- length(seg.sites(mydna))
  L <- length(mydna[1,]) 
  ss_normalized <- S/L
  td <- tajima.test(mydna)
  stats <- data.frame(filename, pi, ss, ss_normalized, 
                      td, td$D, td$Pval.normal, td$Pval.beta)
  df_total <- rbind(df_total,stats)
}

# write dataframe to file
getwd()
write.csv(df_total, "./output/popgenstats.csv")
