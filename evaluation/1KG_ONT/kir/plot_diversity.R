library(ggplot2)


pdf(file="hla/hla_diversity.pdf", width=8, height=4, onefile=FALSE)
df<-read.table("hla/diversity_df.csv", sep=",", header=TRUE)

## maintain the order of the genes
df$Locus <- factor(df$Locus, levels = c('HLA-A', 'HLA-B', 'HLA-C', 'HLA-E', 'HLA-F', 'HLA-G', 'HLA-DRA', 'HLA-DQA1', 'HLA-DQA2', 'HLA-DQB1', 'HLA-DQB2', 'HLA-DPA1', 'HLA-DPA2', 'HLA-DPB1', 'HLA-DPB2', 'HLA-DMA', 'HLA-DMB', 'HLA-DOA', 'HLA-DOB', 'HLA-DRB1', 'HLA-DRB3', 'HLA-DRB4', 'HLA-DRB5', 'HLA-H', 'HLA-J', 'HLA-K', 'HLA-L', 'HLA-N', 'HLA-P', 'HLA-S', 'HLA-T', 'HLA-U', 'HLA-V', 'HLA-W', 'HLA-Y', 'MICA', 'MICB', 'TAP1', 'TAP2'))

# Stacked barplot with multiple groups

ggplot(data=df, aes(x=Locus, y=Diversity)) +
  geom_boxplot(color="#999999") +theme_classic()+
   theme(legend.position = "top", axis.text.x = element_text(size = 8))+
   scale_x_discrete(guide = guide_axis(angle = 90))+
  xlab("")

dev.off()

