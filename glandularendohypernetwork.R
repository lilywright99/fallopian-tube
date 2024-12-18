library(dendextend)
library(gplots)
library(readr)
library(magrittr)
library(tidyverse)
library(ggplot2)
library(ggrepel)

endomatrix<-read_csv('/Users/user/Documents/Fallopiantube/R/glandularendojune.csv')

sd.scores<-apply(endomatrix,2,sd)
endomatrix<-endomatrix[,which(sd.scores>0)]

DEGs<-read_lines ("/Users/user/Documents/Fallopiantube/R/broadDEGsjun.txt",)

endocor<-cor(endomatrix [,na.omit(match(DEGs, colnames(endomatrix)))],
             endomatrix [,-na.omit(match(DEGs, colnames(endomatrix)))])

endobin<-abs(endocor)
endobin[which(endobin>sd(endocor))]<-1
endobin[which(endobin!=1)]<-0 

endohyp<-endobin%*% t(endobin)   

endohm<-heatmap.2(endohyp,trace="none",labRow=F,labCol=F, dendrogram = 'row', col=bluered)

k=6
endodend<-as.hclust(endohm$rowDendrogram)
endodend<-dendextend::colour_branches(endodend,k,groupLabels=T)
plot(endodend)
print(endodend)
endoclusters<-dendextend::cutree(endodend,k=6,order_clusters_as_data=F)


endo_desired_branch <- dendextend::cutree(
  endodend,
  k = 6,
  order_clusters_as_data = FALSE)

table(endo_desired_branch)

branch3<- data.frame(endo_desired_branch[endo_desired_branch == 3 | endo_desired_branch == 4])

endoHNbranches<-data.frame(endo_desired_branch)
write.csv(endoHNbranches,"/Users/user/Documents/Fallopiantube/R/glandularhypernetworkbranches.csv")

#calculating row sums for all DEGs 
rowsum <-rowSums(endobin)

#plot a distribution curve
hist(rowsum)

rowsumdf <- (data.frame(rowsum = rowsum))
rowsumdf

ggplot(rowsumdf, aes(x=rowsum)) + geom_histogram(bins=60)+ geom_density()
ggplot(rowsumdf, aes(x=rowsum)) + geom_density( fill ="lightblue", alpha= 0.5)

gwasgenes<-read_lines('/Users/user/Documents/Fallopiantube/R/broad_gwasgenes.tsv')

rowsum<- subset(endobin, rownames(endobin) %in% gwasgenes)
rowsumDEGS <-rowSums(rowsum)

rowsumDEGS <-data.frame(rowsum = rowsumDEGS, gene = names(rowsumDEGS))

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

write.csv(rowsumrank, '/Users/user/Documents/Fallopiantube/R/endoglandularrowsumrank.csv')

rowsumsd <-sd(rowsumdf$rowsum)
rowsummean <- mean(rowsumdf$rowsum)

rankplot + 
  geom_vline(aes(xintercept = rowsummean), color = "red", linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = rowsummean + rowsumsd), color = "blue", linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = rowsummean - rowsumsd), color = "blue", linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = rowsummean + 2*rowsumsd), color = "blue", linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = rowsummean - 2*rowsumsd), color = "blue", linetype = "dashed", size = 1)
