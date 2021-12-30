collapsed <- read.csv("./collapsed_OTU5.csv", header = TRUE)
collapsed$Strain <- gsub("^.{0,4}", "", collapsed$Strain)

stock <- read.csv("./in-80C.csv", header = TRUE)


merged <- merge(stock, collapsed)

write.csv(merged, "./collapsed_stock.csv")
