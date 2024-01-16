library(Rsamtools) # samtools interface for R
library(GenomicAlignments) # containers for storing gene alignments - read counting, computing coverage, nuc content etc.
library(tidyr) # package for 'tidying' data with limited set of functions for reshaping
library(dplyr) # set of functions for data manipulation
library(edgeR) # differential expression analysis
library(ggplot2)
library(factoextra)
library(ggfortify)
library(viridis)
library(limma)
library(gplots)

getwd()
setwd("/home/talia-backman/Desktop/Ps1524_tailocin")

###############################################
###### for comparing the T1 vs ctrl ##########
###############################################
sample_order_p = read.table("./data/TnSeq/sample_info_p25c2.csv",
                            sep = ",",header = T)
rownames(sample_order_p) <- sample_order_p$sample
sample_order_p_t3 <- sample_order_p %>% filter(time_point=="T1")

#######################
# REMOVE THE BLANKS #
######################
# makes a new table, "counts", by reading in the delimited file (the gffs) - the first row does not contain column headers
counts = read.delim("./data/TnSeq/p25c2_headers_counts.gff", header = T)
# count_table takes all rows and columns 10 on from 'counts' (I understand dim(counts)[2] to give the number of columns in 'counts'). Why is this starting from 10? <- 10 is where the counts start. Before that is information about the region being counted. Is each column after for one of the samples/bams?
count_list <- names(which(colSums(counts[,c(10:dim(counts)[2])])>10000))
count_table = counts[,count_list]

#reorder sample_order
sample_order_p_t3 = sample_order_p_t3[count_list,]
sample_order_p_t3 <- na.omit(sample_order_p_t3)
rownames(count_table) = counts$gene_id

# subset all sample count table to only "t3" samples
count_table_p_t3 <- count_table[, colnames(count_table) %in% sample_order_p_t3$sample]
# take the subsetted t3 count table and 
# creates a DGEList object from a table of counts 
# (rows=features, columns=samples)
dge_p = DGEList(counts=count_table_p_t3, samples = sample_order_p_t3)
dge_p
# duplicate the DGEList for some reason..
dge.trimmed_p_t3 <- dge_p
dge.trimmed_p_t3

M <- median(dge.trimmed_p_t3$samples$lib.size) * 1e-6
#0.318706
L <- mean(dge.trimmed_p_t3$samples$lib.size) * 1e-6
#0.4002404
table(rowSums(dge.trimmed_p_t3$counts==0)==9)
#FALSE  TRUE 
# 5220   114 
lcpm <- cpm(dge.trimmed_p_t3, log=TRUE)
summary(lcpm)
log2(10/M + 2/L)
#5.184831

# create list with mutant library strain and time point
Group <- factor(paste(dge.trimmed_p_t3$samples$mutant_library_strain,dge.trimmed_p_t3$samples$treatment,sep="."))
Group
# make new table with the Group column present
cbind(dge.trimmed_p_t3$samples,Group=Group)

# what is going on here??
keep_1 <- filterByExpr(dge.trimmed_p_t3, Group)
keep_1
dge.trimmed_p__t3_1 <- dge.trimmed_p_t3[keep_1,, keep.lib.sizes=FALSE]
dim(dge.trimmed_p__t3_1)
#[1] 4141   29
keep_2 <- filterByExpr(dge.trimmed_p_t3, Group, min.count = 5)
dge.trimmed_p_t3_2 <- dge.trimmed_p_t3[keep_2,, keep.lib.sizes=FALSE]
dim(dge.trimmed_p_t3_2)
#[1] 4498  29
keep_3 <- filterByExpr(dge.trimmed_p_t3, Group, min.count = 0, min.total.count = 50)
dge.trimmed_p_t3_3 <- dge.trimmed_p_t3[keep_3,, keep.lib.sizes=FALSE]
dim(dge.trimmed_p_t3_3)
#[1] 4799  29
keep_4 <- filterByExpr(dge.trimmed_p_t3, Group, min.count = 0, min.total.count = 50)
dge.trimmed_p_t3_4 <- dge.trimmed_p_t3[keep_4,, keep.lib.sizes=FALSE]
dim(dge.trimmed_p_t3_4)
#[1] 4799  29

dge.trimmed_p_t3_3$samples
dge.trimmed_p_t3_3$samples
dge.trimmed_p <- calcNormFactors(dge.trimmed_p_t3_3)
dge.trimmed_p

samples <- dge.trimmed_p_t3_3$samples
lcpm <- cpm(dge.trimmed_p_t3_3, log=TRUE)
dge_counts <- as.data.frame(lcpm)
dge_counts <- t(dge_counts)
dge.trimmed_p_t3_3$samples

PCA_dge <- prcomp(dge_counts)

plot(PCA_dge)

pdf("./figures/pca.pdf")
autoplot(PCA_dge, data=sample_order_p_t3, colour ='treatment')
dev.off()
PCi<-data.frame(PCA_dge$x,treatment=sample_order_p_t3$treatment)
PCi

pdf("./figures/pca_p25c2_lcpm_t3.pdf", useDingbats = FALSE, fonts = "ArialMT")
ggplot(PCi,aes(x=PC1,y=PC2,col=treatment))+
  geom_point(size=3,alpha=0.5) + #Size and alpha just for fun
  scale_color_manual(values = c("#FDE725FF","#440154FF","#33638DFF"))+ #your colors here
  theme_classic() +
  ylab("PC2 (4.04%)") + xlab("PC1 (88.27%)")
dev.off()

cutoff <- 1
drop <- which(apply(cpm(dge.trimmed_p_t3_3), 1, max) < cutoff)
d <- dge.trimmed_p_t3_3[-drop,]
dim(d)
#0 29
group <- interaction(dge.trimmed_p_t3_3$samples$treatment, dge.trimmed_p_t3_3$samples$time_point)
lcpm <- cpm(dge.trimmed_p_t3_3, log=TRUE)
par(mfrow=c(1,2))
col.group <- group
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
col.group <- 
  as.character(col.group)

pdf("./figures/MDS.pdf")
plotMDS(lcpm, labels=group, col=col.group)
dev.off()


# Specify the model to be fitted. 
# We do this before using voom since voom uses variances of the model residuals 
# (observed - fitted)
design <- model.matrix(~0+Group, data = dge.trimmed_p_t3_3$samples)
# Voom
v <- voom(dge.trimmed_p_t3_3, design, plot = TRUE)
v
# What is voom doing?
#Counts are transformed to log2 counts per million reads (CPM), where “per million reads” is defined based on the normalization factors we calculated earlier
#A linear model is fitted to the log2 CPM for each gene, and the residuals are calculated
#A smoothed curve is fitted to the sqrt(residual standard deviation) by average expression (see red line in plot above)
#The smoothed curve is used to obtain weights for each gene and sample that are passed into limma along with the log2 CPMs.

# fitting linear models in limma
# lmFit fits a linear model using weighted least squares for each gene:
fit <- lmFit(v, design)
head(coef(fit))
# Comparisons between groups (log fold-changes) are obtained as contrasts of these fitted linear models:
# Specify which groups to compare:
# p25c2.WT_tailocin - p25c2.KO_tailocin 
Group
contr_t3 <- makeContrasts(Groupp25c2.WT_tailocin - Groupp25c2.KO_tailocin,
                          levels=colnames(coef(fit)))
contr_t3

#Estimate contrast for each gene
tmp <- contrasts.fit(fit, contr_t3)
# Empirical Bayes smoothing of standard errors (shrinks standard errors that are much larger or smaller than those from other genes towards the average standard error)
tmp <- eBayes(tmp)
summary(decideTests(tmp))
# Groupp25c2.WT_tailocin - Groupp25c2.KO_tailocin
#Down                           4711
#NotSig                           18
#Up                               70

# What genes are most differentially expressed?
genes <- topTable(tmp, sort.by = "P", number = Inf)
genes

# how many "differentially expressed" genes are there?
length(which(genes$adj.P.Val < 0.05)) # 4781 significant genes..

group_p <- factor(paste(sample_order_p_t3$treatment,sample_order_p_t3$time_point,sep="."))
group_p
sample_order_p_t3$group <- group_p

# E values are the equivalent (though the equation differs) to edgeR counts per million. This is like the corrected reads
v_keep <- v$E[rownames(genes),]
v_keep_reorder <- v$E[,sample_order_p_t3[order(sample_order_p_t3$treatment),]$sample]