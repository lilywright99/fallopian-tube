library(dendextend)
library(gplots)
library(readr)
library(magrittr)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ComplexHeatmap)

FT<-read_csv('/Users/user/Documents/Fallopiantube/R/secretoryFTjune.csv')

sd.scores<-apply(FT,2,sd)
FT<-FT[,which(sd.scores>0)]

DEGs<-read_lines ("/Users/user/Documents/Fallopiantube/R/broadDEGsjun.txt",)

cor_data<-cor(FT [,na.omit(match(DEGs, colnames(FT)))],
                 FT [,-na.omit(match(DEGs, colnames(FT)))])

bin<-abs(cor_data)
bin[which(bin>sd(cor_data))]<-1
bin[which(bin!=1)]<-0 

hyp<-bin%*% t(bin)   

hm<-heatmap.2(hyp,trace="none",labRow=F,labCol=F, dendrogram = 'row', col='bluered')

k=4
dend<-as.hclust(hmv3$rowDendrogram)
dend<-dendextend::colour_branches(dend,k,groupLabels=T)
plot(dend)
print(dend)
clusters<-dendextend::cutree(dend,k=4,order_clusters_as_data=F)

desired_branch <- dendextend::cutree(
  dend,
  k = 4,
  order_clusters_as_data = FALSE)

table(desired_branch)

branch78<-data.frame(desired_branch[desired_branch == 1 | desired_branch == 2])


desired_branch<-data.frame(desired_branch)
write.csv(desired_branch,"/Users/user/Documents/Fallopiantube/R/FThypernetworkbranches.csv")

#calculating row sums for all DEGs 
rowsum <-rowSums(bin)

#plot a distribution curve
hist(rowsum)

rowsumdf <- (data.frame(rowsum = rowsum))
rowsumdf
write.csv(rowsumdf,"FTDEGrowsums.csv")

ggplot(rowsumdf, aes(x=rowsum)) + geom_histogram(bins=60)+ geom_density()
ggplot(rowsumdf, aes(x=rowsum)) + geom_density( fill ="lightblue", alpha= 0.5)

gwasgenes<-read_lines('/Users/user/Documents/Fallopiantube/R/broad_gwasgenes.tsv')

rowsum<- subset(bin, rownames(bin) %in% gwasgenes)
rowsumDEGS<-rowSums(rowsum)

rowsumDEGS <- data.frame(rowsum = rowsumDEGS, gene = names(rowsumDEGS))

#calculate row sums for gwas genes in FT 
density_plot <- ggplot(data = data.frame(rowsum = rowsumdf), aes(x = rowsum)) +
  geom_density(fill = "lightblue", alpha = 0.5) +
  theme_minimal() + theme( panel.grid.major = element_blank(), 
                           panel.grid.minor = element_blank()) + labs(title = "Density Plot of DEG Row Sums", x = "Row Sum", y = "Density")

density_plot

HNDEGSplot <- density_plot + geom_vline(data = rowsumDEGS, aes(xintercept = rowsum), color = "darkblue", linetype = "dashed")

HNDEGSplotgenes <- HNDEGSplot + geom_text_repel(data = rowsumDEGS, aes(x = rowsum, y = 0.0010, label = gene), 
                                                angle = 90, vjust = 1.5, hjust = 0, size = 3, color = "black")

HNDEGSplotgenes

#annotate with rank

rowsumdf$Rank <- rank(rowsumdf$rowsum, ties.method = "average")
rowsumdf <- rowsumdf[order(rowsumdf$Rank), ]
rowsumdf <- data.frame(rowsumdf, gene = rownames(rowsumdf))

#crop to ranks of gwas genes only to plot
rowsumrank <-  subset(rowsumdf, rownames(rowsumdf) %in% gwasgenes)

rankplot <- HNDEGSplot +
  geom_text_repel(data = rowsumrank, 
                  aes(x = rowsum, y = 0.00075, label = paste(gene, "\n", 'R=', Rank)), 
                  angle = 90, vjust = 0.5, hjust = 1, size = 3, color = "black", 
                  nudge_y = 0.0002, force = 1, direction = "y")
rankplot

write.csv(rowsumrank, '/Users/user/Documents/Fallopiantube/R/FTrowsumrank.csv')

rowsumsd <-sd(rowsumdf$rowsum)
rowsummean <- mean(rowsumdf$rowsum)

rankplot + 
  geom_vline(aes(xintercept = rowsummean), color = "red", linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = rowsummean + rowsumsd), color = "blue", linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = rowsummean - rowsumsd), color = "blue", linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = rowsummean + 2*rowsumsd), color = "blue", linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = rowsummean - 2*rowsumsd), color = "blue", linetype = "dashed", size = 1)

