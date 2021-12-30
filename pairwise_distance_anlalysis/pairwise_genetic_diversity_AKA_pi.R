# this workflow was derived from: https://evolutionarygenetics.github.io/Chapter7.html

# load libraries
library(ape)
library(PopGenome)
library(pegas)
library(tidyverse)

# FOR BASEPLATE J-LIKE PROTEIN
#Reading sequence data into R.
mydna <- read.dna("./dat/baseplate_jlike_protein.pal2nal", format = "fasta")

#Exploring DNA sequence data
myalign <- as.alignment(mydna)
myalign$seq
alview(mydna)

####################################################################
###########Calculating basic sequence statistics ###################
####################################################################

#Base composition
base.freq(mydna) #a= 18% c=32% g=33% t=17%
GC.content(mydna) # 65%
sum(base.freq(mydna)[c(2, 3)]) # 65%

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
nuc.div(mydna) # 0.03826307, already standardized

#Sample size and sequence statistics
#Segregating sites
# use sapply to loop
ss <- sapply(2:1390, function(z){
  length(seg.sites(mydna[1:z, ]))
})
# plot figure
plot(2:1390, ss, col = "red", xlab = "No of sequences", ylab = "Segregating sites", las = 1)
#the number of segregating sites is biased by the number of sequences we include in the data. It increases our probability of observing 
#a polymorphism and also since all polymorphisms are given equal weighting in the calculation of the number of segregating sites, any 
#polymorphism will increase the value by 1.
#Nucleotide Diversity (pi)
# use sapply to loop
nd <- sapply(2:1390, function(z){
  nuc.div(mydna[1:z, ])
})
plot(2:1390, nd, col = "blue", xlab = "No of sequences", ylab = expression(pi), las = 1)


# Inferring evolutionary processes using Tajima’s D
tajima.test(mydna) # Tajima's D = -1.477048, p=values = 0.1396628, 0.106838








