library(dendextend)
library(gplots)
library(readr)
library(magrittr)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ggdendro)
library(reshape2)

gwasstats<- read.csv("~/Documents/Fallopiantube/R/3NHSbroad_gwasranks.csv", header = TRUE)
rankedgenes <- read.table("~/Documents/Fallopiantube/R/rowsumranksplot.csv", sep = ",", header=TRUE)

gwasstats<- subset(gwasstats, gwasstats$gene %in% c(rankedgenes$fallopian.tube.gene, rankedgenes$luminal.gene, rankedgenes$glandular.gene))

ggplot(gwasstats, aes(x = factor(gene, level=c('PGGHG','ANKRD36B','ANKRD36','SLC7A2','MUC1','CLDN1','HLA-DRB1','HLA-DRB5','PEX6','PLXNA4','NR2F1','GLS','CLGN')), y = Rank)) +
  geom_point(aes(color = Hypernetwork), position = position_dodge(width = 0), size = 1.5) + 
  theme_minimal() +
  labs(
    title = "Ranks of Ectopic Pregnancy Associated Genes
    in Fallopian tube and Endometrial Hypernetworks",
    x = "Gene",
    y = "Rank",
    color = "Hypernetwork") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 14, hjust = 0.5))
