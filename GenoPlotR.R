# this script doesn't work, keep in case of future use
# load libraries
library(genoPlotR)

# read in data for subcluster 6 within cluster 2
p25.A12 <- read_dna_seg_from_file("./R_dat/ffns/136.plate25.A12.pilon.contigs_renamed.phages_lysogenic.ffn")
p26.D6 <- read_dna_seg_from_file("./R_dat/ffns/403.plate26.D6.pilon.contigs_renamed.phages_lysogenic.ffn")
p26.F9 <- read_dna_seg_from_file("./R_dat/ffns/92_.plate26.F9.pilon.contigs_renamed.phages_lysogenic.ffn")
p25.D2 <- read_dna_seg_from_file("./R_dat/ffns/662.plate25.D2.pilon.contigs_renamed.phages_lysogenic.ffn")

read_dna_seg_from_genbank("./R_dat/gbks/plate25.A12.pilon.contigs_renamed.phages_combined (copy).gbk")


p25.A12gbk <- read_dna_seg_from_genbank(file = "./R_dat/gbks/136_plate25.A12.pilon.contigs_renamed.phages_combined.temp.gbk", tagsToParse = "CDS")

read_dna_seg_from_file("./R_dat/136_p25.A12.gbk")



#Import your  Comparison files.
p25.A12_p26.D6.comparison <- read_comparison_from_blast("./R_dat/fnas/p25.A12_p26.D6.blastn",
                                                        sort_by = "per_id",filt_high_evalue = NULL,
                                                        filt_low_per_id = NULL, filt_length = NULL,
                                                        color_scheme = NULL)
p25.A12_p26.F9.comparison <- read_comparison_from_blast("./R_dat/fnas/p25.A12_p26.F9.blastn",
                                                        sort_by = "per_id",filt_high_evalue = NULL,
                                                        filt_low_per_id = NULL, filt_length = NULL,
                                                        color_scheme = NULL)
p25.A12_p25.D2.comparison <- read_comparison_from_blast("./R_dat/fnas/p25.A12_p25.D2.blastn",
                                                        sort_by = "per_id",filt_high_evalue = NULL,
                                                        filt_low_per_id = NULL, filt_length = NULL,
                                                        color_scheme = NULL)


#check that your comparison file is acceptable:
?is.comparison
is.comparison(p25.A12_p26.D6.comparison) # TRUE, Yay!!
is.comparison(p25.A12_p26.F9.comparison) # TRUE
is.comparison(p25.A12_p25.D2.comparison) # TRUE

#Give you dna_seg sequences CDS annotation
p25.A12_annot <-auto_annotate(p25.A12)
p26.D6_annot <-auto_annotate(p26.D6)
p26.F9_annot <- auto_annotate(p26.F9)
p25.D2_annot <- auto_annotate(p25.D2)

is.annotation(p25.A12_annot) # TRUE

?annotation()

# Make your DNA-comparison image

xlims1 <- list(c(2,1615),
               c(1746,2042),
               c(2039,2386),
               c(2454,3950),
               c(3934,4155),
               c(4152,4742),
               c(4829,5167),
               c(5148,5537),
               c(5842,6279),
               c(6422,7066),
               c(7215,7472),
               c(7561,7743),
               c(7790,8032))


plot_gene_map(dna_segs=list(p25.A12, p26.D6, p26.F9, p25.D2), comparisons = list(p25.A12_p26.D6.comparison, p25.A12_p26.F9.comparison, p25.A12_p25.D2.comparison), 
              annotations=list(p25.A12_annot, p26.D6_annot, p26.F9_annot, p25.D2_annot),
              dna_seg_labels = c("p25.A12","p26.D6","p26.F9","p25.D2"),main = "Subcluster 6 gene annotation")



plot_gene_map(dna_segs=list(p25.A12, p26.D6, p26.F9, p25.D2), comparisons = list(p25.A12_p26.D6.comparison, p25.A12_p26.F9.comparison, p25.A12_p25.D2.comparison),
              dna_seg_labels = c("p25.A12","p26.D6","p26.F9","p25.D2"),main = "Subcluster 6 gene annotation")


?plot_gene_map
