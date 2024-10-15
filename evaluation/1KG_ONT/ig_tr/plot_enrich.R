library(ggplot2)
library(dplyr)
require(gridExtra)
library(grid)
library(cowplot)
library(ggpubr)
# library(ggstance)

df<-read.table("vdj_compare_hete_num_AFR.csv", sep=",", header=TRUE)


df <- filter(df, p.adj < 0.05)

pdf(file="vdj_comp_hete_enriched_AFR.pdf", width=6, height=10, onefile=FALSE)


# Create the plot
h1<-ggplot(df, aes(x = Gene, y = fold_change)) +
  geom_point(aes(size = total, color =  log10(p.adj))) +
  scale_color_gradient(low ="#2a347a", high = "#d6d69b") +
  theme_bw()+
  coord_flip()+
  xlab("Gene")+
  ylab("Fold Change")+
  scale_y_continuous(trans='log2')+
  geom_hline(yintercept = 1, linetype = "dashed")
h1

dev.off()