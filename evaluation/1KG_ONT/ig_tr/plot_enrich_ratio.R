library(ggplot2)
library(dplyr)
require(gridExtra)
library(grid)
library(cowplot)
library(ggpubr)
# library(ggstance)

df<-read.table("vdj_compare_Enriched.csv", sep=",", header=TRUE)



pdf(file="vdj_compare_Enriched.pdf", width=4, height=4, onefile=FALSE)

p<-ggplot(df, aes(x=Pop, y=loci_num, fill=Type)) +
  geom_bar(stat="identity")+theme_minimal()+
  theme_classic()+
  ylab("No. of loci")+xlab("Population")
  # geom_text(aes(label=loci_num), vjust=1.6, color="white",
  #           position = position_dodge(0.9), size=3.5)

p+scale_fill_brewer(palette="Dark2")


dev.off()