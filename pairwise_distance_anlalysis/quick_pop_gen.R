# load libraries
library(ape)
library(PopGenome)
library(pegas)
library(tidyverse)

mydna <- read.dna("./otu5_lowdivergentcluster_dat/VOG01867.pal2nal", format = "fasta")

base.freq(mydna) #a= 18% c=32% g=33% t=17%
GC.content(mydna) # 65%

#Segregating sites
seg.sites(mydna)
length(seg.sites(mydna)) #635 positions where polymorphisms occur, but this is not standardised to the length of the sequences
# to achieve standardisation we do the following:
# get segregating sites
S <- length(seg.sites(mydna))
# set sequence length
L <- 1035
# standardise segregating sites by sequence length
s <- S/L 
s # 0.6135266

#Nucleotide diversity, or the average number of differences between sequences in a population or sample.
nuc.div(mydna)

# Inferring evolutionary processes using Tajima’s D
tajima.test(mydna)
