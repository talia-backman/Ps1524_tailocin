# load libraries
library(factoextra)
# set the seed
set.seed(123)

# load and tidy phage data from mash 
tbl <- read.table("./data/mash_data/final_mash_results.txt",header=FALSE,row.names=1)
colnames(tbl) <- rownames(tbl)

# find optimal amount of viral element clusters in our mash data
# Elbow method
fviz_nbclust(tbl, kmeans, method = "wss",k.max = 10) +
 geom_vline(xintercept = 4, linetype = 2) +
  labs(subtitle = "Elbow method") # k = 4
# Silhouette method
fviz_nbclust(tbl, kmeans, method = "silhouette",k.max=10)+
 labs(subtitle = "Silhouette method") # k = 2

# plot k-means using 4 clusters 
tbl2 <- as.matrix(tbl)
p4 <- kmeans(tbl2, 4,nstart=20)
p2 <- fviz_cluster(p4, data = tbl2, ellipse = FALSE,shape = "circle",palette = c("#58137b", "#1c78ac","#1ccf87","#ece715"), geom="point", 
                   ggtheme = theme_classic(),
                   main="") + theme(legend.position = "none", 
                                    axis.text = element_text(size=20),
                                    axis.title = element_text(size=24))
# save clustering figure
pdf("./figures/01_phage_annotation/clustering.pdf")
p2
dev.off()