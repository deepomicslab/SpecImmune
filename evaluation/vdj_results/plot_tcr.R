
library(ggplot2)
library(cowplot)
# install.packages("gridExtra")
# require(gridExtra)
# library(grid)


pdf(file="figures/tcr11.pdf", width=12, height=8, onefile=FALSE)


df<-read.table("tcr_11.csv", sep=",", header=TRUE)


p1<-ggplot(data=df, aes(x=gene, y=accuracy)) +
  geom_bar(stat="identity", width=0.5, fill="#E69F00")+
  theme_minimal()+
   theme(legend.position = "none", axis.text.x = element_text(size = 6))+
   scale_x_discrete(guide = guide_axis(angle = 90))+
  xlab("")+
  ylab("Accuracy")
# p1

p2<-ggplot(data=df, aes(x=gene, y=total)) +
  geom_bar(stat="identity", width=0.5, fill="#E69F00")+
  theme_minimal()+
   theme(legend.position = "none", axis.text.x = element_text(size = 6))+
   scale_x_discrete(guide = guide_axis(angle = 90))+
  xlab("TCR Loci")+
  ylab("No. of alleles")
# p1

prow <- plot_grid(
  p1 + theme(legend.position="none"),
  p2 + theme(legend.position="none"),
  align = 'vh',
  hjust = -1,
  nrow = 2
)

plot_grid(prow, legend, ncol = 1, rel_heights = c(1, .1))


# df <- arrange(df, desc(fold))
# my_order = c("SCFA", "Butyrate", "Butyryl-CoA", "propionate", "Succinate", "Pyruvate", "Acetyl-CoA", "Lactate", "Acetate")
# df$scfa <- factor(df$scfa, levels = my_order)
#df <- filter(df, p.adj < 0.05)
# ggplot(df,aes(y=accuracy,x=gene)+ 
#   geom_bar(stat="identity") +
#   theme_bw()+
# #   coord_flip()+
#   xlab("")+
#   ylab("Fold Enrichment"))#+
  #scale_y_continuous(trans='log2')+
  #geom_hline(yintercept = 1, linetype = "dashed") +
#   theme(strip.text.y = element_text(angle = 0), legend.position = "bottom", axis.text.y = element_text(size = 15)) +
#   geom_text(aes(label=paste("p =", format(p.adj, scientific = TRUE, digits = 3))), 
#             position=position_dodge(width=0.9), 
#             vjust=0.5,hjust=1.1, size=3, color="white")



dev.off()