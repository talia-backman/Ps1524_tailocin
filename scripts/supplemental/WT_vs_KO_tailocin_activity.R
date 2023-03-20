# load libraries
library(ggplot2)
library(tidyr)
library(dplyr)
library(latex2exp)
library(cowplot)
library(ggpubr)

# read in data
# WT vs KO tecan data (OD600 over time)
dat <- read.csv("./data/supplemental_data/WT_vs_KO_tailocin_activity.csv",header=T,sep=",")

# clean dat
# transform dataframe from wide to long format
df.long <- pivot_longer(dat, cols=2:5, names_to = "Sample", values_to = "OD600")
# divide seconds by 60 to get minutes, then hours
df.long$Time_min <- df.long$Time..s./60
df.long$Time_hours <- df.long$Time_min/60
df.long <- na.omit(df.long)
# Remove media only from the data
df.long <- subset(df.long, Sample != "Media")

# plot data
p <- ggplot(df.long, aes(x=Time_hours, y=OD600, color=Sample)) + geom_point() + 
  geom_smooth(aes(group=Sample),se=FALSE) + 
  theme_bw() + xlab("Time after treatment (Hours)") + ylab(bquote(OD[600])) +
  scale_color_manual(values=c("#58137b", "#1c78ac","#1ccf87","#ece715")) 

# export as pdf
pdf("./figures/supplemental/WT_vs_KO_tailocin_activity.pdf")
p
dev.off()
