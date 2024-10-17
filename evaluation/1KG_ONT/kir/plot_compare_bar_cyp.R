library(ggplot2)
library(ggpubr)
library(ggplot2)
library(grid)
library(cowplot)

pdf(file="figures/hla_class_kir_cyp.pdf", width=4, height=4, onefile=FALSE)
df<-read.table("tree/cyp_hla_kir.csv", sep=",", header=TRUE)

## maintain the order of the genes
# df$Gene <- factor(df$Gene, levels = c('HLA-A', 'HLA-B', 'HLA-C', 'HLA-E', 'HLA-F', 'HLA-G', 'HLA-DRA', 'HLA-DQA1', 'HLA-DQA2', 'HLA-DQB1', 'HLA-DQB2', 'HLA-DPA1', 'HLA-DPA2', 'HLA-DPB1', 'HLA-DPB2', 'HLA-DMA', 'HLA-DMB', 'HLA-DOA', 'HLA-DOB', 'HLA-DRB1', 'HLA-DRB3', 'HLA-DRB4', 'HLA-DRB5', 'HLA-H', 'HLA-J', 'HLA-K', 'HLA-L', 'HLA-N', 'HLA-P', 'HLA-S', 'HLA-T', 'HLA-U', 'HLA-V', 'HLA-W', 'HLA-Y', 'MICA', 'MICB', 'TAP1', 'TAP2'))

# Stacked barplot with multiple groups

# ### sort the df by ALD in reverse order
df <- df[order(df$ALD, decreasing = TRUE),]
df$CYP <- factor(df$CYP, levels = df$CYP)

### select the df with gene=="HLA-DPB1"


vysg<-ggplot(data=df, aes(x=CYP, y=ALD)) +
  geom_bar(stat="identity",fill="#9fc3d5")+theme_bw()+
   theme(legend.position = "top", axis.text.x = element_text(size = 8))+
   scale_x_discrete(guide = guide_axis(angle = 90))+
  xlab("")

vysg<-vysg+facet_wrap(~ type , ncol = 1, scales = "free_x")+
xlab('')+
ylab('minALD')+
  theme(legend.position = "none")

vysg










# df2 <- df[df$type == "HLA",]
# df2 <- df2[order(df2$ALD, decreasing = TRUE),]
# df2$CYP <- factor(df2$CYP, levels = df2$CYP)

# df3 <- df[df$type == "KIR",]
# df3 <- df3[order(df3$ALD, decreasing = TRUE),]
# df3$CYP <- factor(df3$CYP, levels = df3$CYP)

# p1<-ggplot(data=df2, aes(x=CYP, y=ALD)) +
#   geom_bar(stat="identity",fill="#9fc3d5")+theme_bw()+
#    theme(legend.position = "top", axis.text.x = element_text(size = 8))+
#    scale_x_discrete(guide = guide_axis(angle = 90))+
#   xlab("")

# p2<-ggplot(data=df3, aes(x=CYP, y=ALD)) +
#   geom_bar(stat="identity",fill="#9fc3d5")+theme_bw()+
#    theme(legend.position = "top", axis.text.x = element_text(size = 8))+
#    scale_x_discrete(guide = guide_axis(angle = 90))+
#   xlab("")


# prow <- plot_grid(
#   p1 + theme(legend.position="none"),
#   p2 + theme(legend.position="none"),
#   align = 'vh',
#   hjust = -1,
#   nrow = 2
# )
# legend <- get_legend(
#   p1 +
#     guides(color = guide_legend(nrow = 1)) +
#     theme(legend.position = "bottom")
# )
# # legend = get_plot_component(p4, 'guide-box-top', return_all = TRUE)
# plot_grid(prow, legend, ncol = 1, rel_heights = c(1, .1))

# dev.off()