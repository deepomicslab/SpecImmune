library(ggplot2)
library(ggpubr)

pdf(file="figures/hla_class_kir2.pdf", width=6, height=3, onefile=FALSE)
df<-read.table("tree/hla_class_kir.csv", sep=",", header=TRUE)

## maintain the order of the genes
# df$Gene <- factor(df$Gene, levels = c('HLA-A', 'HLA-B', 'HLA-C', 'HLA-E', 'HLA-F', 'HLA-G', 'HLA-DRA', 'HLA-DQA1', 'HLA-DQA2', 'HLA-DQB1', 'HLA-DQB2', 'HLA-DPA1', 'HLA-DPA2', 'HLA-DPB1', 'HLA-DPB2', 'HLA-DMA', 'HLA-DMB', 'HLA-DOA', 'HLA-DOB', 'HLA-DRB1', 'HLA-DRB3', 'HLA-DRB4', 'HLA-DRB5', 'HLA-H', 'HLA-J', 'HLA-K', 'HLA-L', 'HLA-N', 'HLA-P', 'HLA-S', 'HLA-T', 'HLA-U', 'HLA-V', 'HLA-W', 'HLA-Y', 'MICA', 'MICB', 'TAP1', 'TAP2'))

# Stacked barplot with multiple groups

### sort the df by ALD in reverse order
df <- df[order(df$ALD, decreasing = TRUE),]

### select the df with gene=="HLA-DPB1"
df2 <- df[df$KIR == "KIR2DL1",]
### maintain the order
df$gene <- factor(df$gene, levels = df2$gene)

vysg<-ggplot(data=df, aes(x=gene, y=ALD)) +
  geom_bar(stat="identity",fill="#9fc3d5")+theme_bw()+
   theme(legend.position = "top", axis.text.x = element_text(size = 8))+
   scale_x_discrete(guide = guide_axis(angle = 90))+
  xlab("")

vysg<-vysg+facet_wrap( ~ KIR , ncol = 4)+
xlab('')+
ylab('minALD')+
  theme(legend.position = "none")

vysg

dev.off()