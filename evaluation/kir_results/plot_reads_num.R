library(ggplot2)
library(grid)
library(cowplot)
library(dplyr)


pdf(file="figures/read_num.pdf", width=16, height=6, onefile=FALSE)

df<-read.table("HGSCV2_hifi_read_num.csv", sep=",", header=TRUE)


# ## add a column to the data frame
# df$dataset <- "HPRC"
# df2$dataset <- "HGSCV2"
# ## combine the two data frames
# df <- rbind(df, df2)

p1<- ggplot(data=df, aes(x=gene, y=read_num)) +
  geom_boxplot(color="#E69F00")+
  theme_classic()+
   theme(legend.position = "none", axis.text.x = element_text(size = 10))+
   scale_x_discrete(guide = guide_axis(angle = 90))+
  xlab("")+
  ggtitle('HPRC HiFi')+
  ylab("Count")+scale_color_brewer(palette="Dark2")

p1+ geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)
# prow <- plot_grid(
#   p1 + theme(legend.position="none"),

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
dev.off()