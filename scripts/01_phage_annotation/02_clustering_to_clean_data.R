# load libraries
library(ape)
library(dplyr)

# read in and tidy data
# phage data from mash 
tbl <- read.table("./data/mash_data/final_mash_results.txt",header=FALSE,row.names=1)
colnames(tbl) <- rownames(tbl)
tbl2 <- as.matrix(tbl)
p4 <- kmeans(tbl2, 4,nstart=20)

# clean the cluster data to be able to map them to the pseudomonas phylogeny
df <- data.frame(p4$cluster)
df$genome_ID <- row.names(df)
df[] <- lapply(df, gsub, pattern='late', replacement='')
df$genome_ID <- gsub("^.{0,4}", "", df$genome_ID)
df[] <- lapply(df, gsub, pattern='.fna', replacement='')
rownames(df) <- NULL
df <- df[,c(2,1)]

##Add a column to clusters: Phage #
df['phages'] <- NA
df <- df %>% group_by(genome_ID) %>% dplyr::mutate(phages = row_number())

# subset df based off phages column, rename their cluster column
df1 <- subset(df, phages == 1) #1582
names(df1)[2] <- "phage1"
df2 <- subset(df, phages == 2) # 884
names(df2)[2] <- "phage2"
df3 <- subset(df, phages == 3) # 443
names(df3)[2] <- "phage3"
df4 <- subset(df, phages == 4) # 149
names(df4)[2] <- "phage4"
df5 <- subset(df, phages == 5) # 40
names(df5)[2] <- "phage5"
df6 <- subset(df, phages == 6) # 4
names(df6)[2] <- "phage6"
df7 <- subset(df, phages == 7) # 1
names(df7)[2] <- "phage7"
df8 <- subset(df, phages == 8) # 1
names(df8)[2] <- "phage8"

# make comprehensive phage clustering dataframe with the corresponding genome as rownames
# make all 15 dataframes 1510 observations long but with NA in places that don't have values
temp <- merge(df1, df2, by="genome_ID", all=TRUE)
temp <- merge(temp,df3, by="genome_ID", all=TRUE)
temp <- merge(temp,df4, by="genome_ID", all=TRUE)
temp <- merge(temp,df5, by="genome_ID", all=TRUE)
temp <- merge(temp,df6, by="genome_ID", all=TRUE)
temp <- merge(temp,df7, by="genome_ID", all=TRUE)
all_df <- merge(temp,df8, by="genome_ID", all=TRUE)
all_df <- all_df[ ,c(1,2,4,6,8,10,12,14,16)]

# convert to matrix 
all_df_matrix <- as.matrix(all_df)
rownames(all_df_matrix) <- all_df_matrix[,1]
all_df_matrix <- all_df_matrix[ ,c(2:9)]

# write cleaned dataframe to csv
write.csv(all_df_matrix, "./data/mash_data/clustering_dat_clean.csv")
# manually sort them in excel (ie. make VC1 in column 1 when present, then VC2 in column 2 when present)