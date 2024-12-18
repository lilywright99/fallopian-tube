library(dplyr)
library(ggplot2)
library(ggrepel)
library(tidyverse)

ORA<-read.table("~/Documents/Fallopiantube/R/ORAresults.csv", sep=",", header = TRUE)

ggplot(ORA, aes(x = Tissue, y = Description)) +
  geom_point(aes(size= -log10(FDR), color = -log10(P.Value))) +
  scale_color_gradient(low = "blue", high = "red", name = "-log10(pValue)") +
  scale_size_continuous(range = c(0.2, 15)) +
  theme_minimal() +
  labs(
    title = "GeneOntology Over Representation Analysis 
    for Fallopian Tube and Endometrial Secretory Epithelial Cells",
    x = "Cell Types",
    y = "Pathway Description" ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
   plot.title = element_text(size = 14, face = "bold", hjust = 0.8))





