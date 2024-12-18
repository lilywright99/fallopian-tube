library(dendextend)
library(gplots)
library(readr)
library(magrittr)
library(tidyverse)
library(ggplot2)
library(ggrepel)

lumenalmatrix<-read_csv('/Users/user/Documents/Fallopiantube/R/lumenalendojune.csv')

sd.scores<-apply(lumenalmatrix,2,sd)
lumenalmatrix<-lumenalmatrix[,which(sd.scores>0)]

DEGs<-read_lines ("/Users/user/Documents/Fallopiantube/R/broadDEGsjun.txt",)

lumenalcordata<-cor(lumenalmatrix [,na.omit(match(DEGs, colnames(lumenalmatrix)))],
                    lumenalmatrix [,-na.omit(match(DEGs, colnames(lumenalmatrix)))])

lumenalbin<-abs(lumenalcordata)
lumenalbin[which(lumenalbin>sd(lumenalcordata))]<-1
lumenalbin[which(lumenalbin!=1)]<-0 

lumenalhyp<-lumenalbin%*% t(lumenalbin)   

lumenalhm<-heatmap.2(lumenalhyp,trace="none",labRow=F,labCol=F, dendrogram = 'row', col='bluered')

upper_tri_lumenalhyp <- lumenalhyp
upper_tri_lumenalhyp[row(upper_tri_lumenalhyp) >= col(upper_tri_lumenalhyp)] <- NA


k=7
lumenaldend<-as.hclust(lumenalhm$rowDendrogram)
lumenaldend<-dendextend::colour_branches(lumenaldend,k,groupLabels=T)
plot(lumenaldend)
print(lumenaldend)
lumenalclusters<-dendextend::cutree(lumenaldend,k=9,order_clusters_as_data=F)


lumenal_desired_branch <- dendextend::cutree(
  lumenaldend,
  k = 7,
  order_clusters_as_data = FALSE)

table(lumenal_desired_branch)

branch1<-data.frame(lumenal_desired_branch[lumenal_desired_branch == 4 | lumenal_desired_branch == 6| lumenal_desired_branch == 7])

clustergenes<- list(rownames(branch1))
clustergenes
lumenalHNbranches<-data.frame(lumenal_desired_branch)
write.csv(lumenalHNbranches,"/Users/user/Documents/Fallopiantube/R/lumenalhypernetworkbranches.csv")

#calculating row sums for all DEGs 
rowsum <-rowSums(lumenalbin)

#plot a distribution curve
hist(rowsum)

#create df of rowrums of DEGs 
rowsumdf <- (data.frame(rowsum = rowsum))
rowsumdf

ggplot(rowsumdf, aes(x=rowsum)) + geom_histogram(bins=60)+ geom_density()
ggplot(rowsumdf, aes(x=rowsum)) + geom_density( fill ="lightblue", alpha= 0.5)

gwasgenes<-read_lines('/Users/user/Documents/Fallopiantube/R/broad_gwasgenes.tsv')

#calculate row sums of GWAS genes 
rowsum<- subset(lumenalbin, rownames(lumenalbin) %in% gwasgenes)
rowsumDEGS<-rowSums(rowsum)

rowsumDEGS <- data.frame(rowsum = rowsumDEGS, gene = names(rowsumDEGS))

#calculate row sums for gwas genes in FT 
density_plot <- ggplot(data = data.frame(rowsum = rowsumdf), aes(x = rowsum)) +
  geom_density(fill = "lightblue", alpha = 0.5) +
  theme_minimal() + theme( panel.grid.major = element_blank(), 
                           panel.grid.minor = element_blank()) + labs(title = "Density Plot of DEG Row Sums", x = "Row Sum", y = "Density")

density_plot

HNDEGSplot <- density_plot + geom_vline(data = rowsumDEGS, aes(xintercept = rowsum), color = "darkblue", linetype = "dashed")

HNDEGSgenesplot <- HNDEGSplot + geom_text_repel(data = rowsumDEGS, aes(x = rowsum, y = 0.0015, label = gene), 
                                                angle = 90, vjust = 1.5, hjust = 0, size = 3, color = "black")

HNDEGSgenesplot

#annotate with rank

rowsumdf$Rank <- rank(rowsumdf$rowsum, ties.method = "average")
rowsumdf <- rowsumdf[order(rowsumdf$Rank), ]
rowsumdf <- data.frame(rowsumdf, gene = rownames(rowsumdf))

#crop to ranks of gwas genes only to plot
rowsumrank <-  subset(rowsumdf, rownames(rowsumdf) %in% gwasgenes)

rankplot <- HNDEGSplot +
  geom_text_repel(data = rowsumrank, 
                  aes(x = rowsum, y = 0.002, label = paste(gene, "\n", 'R=', Rank)), 
                  angle = 90, vjust = 0.5, hjust = 1, size = 3, color = "black", 
                  nudge_y = 0.0002, force = 1, direction = "y")
rankplot

write.csv(rowsumrank, '/Users/user/Documents/Fallopiantube/R/endolumenalrowsumrank.csv')

rowsumsd <-sd(rowsumdf$rowsum)
rowsummean <- mean(rowsumdf$rowsum)

rankplot + 
  geom_vline(aes(xintercept = rowsummean), color = "red", linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = rowsummean + rowsumsd), color = "blue", linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = rowsummean - rowsumsd), color = "blue", linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = rowsummean + 2*rowsumsd), color = "blue", linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = rowsummean - 2*rowsumsd), color = "blue", linetype = "dashed", size = 1)

