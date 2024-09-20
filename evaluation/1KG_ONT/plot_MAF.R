library(ggplot2)


pdf(file="hla/hla_MAF.pdf", width=8, height=4, onefile=FALSE)
df<-read.table("hla/hla_maf.csv", sep=",", header=TRUE)

## maintain the order of the genes
df$Gene <- factor(df$Gene, levels = c('HLA-A', 'HLA-B', 'HLA-C', 'HLA-E', 'HLA-F', 'HLA-G', 'HLA-DRA', 'HLA-DQA1', 'HLA-DQA2', 'HLA-DQB1', 'HLA-DQB2', 'HLA-DPA1', 'HLA-DPA2', 'HLA-DPB1', 'HLA-DPB2', 'HLA-DMA', 'HLA-DMB', 'HLA-DOA', 'HLA-DOB', 'HLA-DRB1', 'HLA-DRB3', 'HLA-DRB4', 'HLA-DRB5', 'HLA-H', 'HLA-J', 'HLA-K', 'HLA-L', 'HLA-N', 'HLA-P', 'HLA-S', 'HLA-T', 'HLA-U', 'HLA-V', 'HLA-W', 'HLA-Y', 'MICA', 'MICB', 'TAP1', 'TAP2'))

# Stacked barplot with multiple groups

ggplot(data=df, aes(x=Gene, y=Frequency, fill=Group)) +
  geom_bar(stat="identity")+scale_fill_brewer(palette="Dark2") +theme_classic()+
   theme(legend.position = "top", axis.text.x = element_text(size = 8))+
   scale_x_discrete(guide = guide_axis(angle = 90))+
  xlab("")

dev.off()


## extract the elements with Group as rare
df1 <- df[df$Group == "rare",]
## sort df1 by Frequency
df1 <- df1[order(df1$Frequency),]
df1