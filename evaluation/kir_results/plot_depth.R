library(ggplot2)
library(grid)
library(cowplot)
library(dplyr)


pdf(file="figures/kir_depth_change.pdf", width=8, height=3, onefile=FALSE)

df<-read.table("all_loci_depth.csv", sep=",", header=TRUE)
df1 <- dplyr::filter(df, field == 1)


p1<- ggplot(data=df1, aes(x=depth, y=accuracy, color=dataset)) +
  geom_bar(stat="identity", position=position_dodge())+ #"#E69F00"
  # theme_minimal()+
  theme_classic()+
  geom_line()+
  geom_point()+
  xlab("Depth")+
  ylim(c(0.8,1))+
  ggtitle('1-field')+
  ylab("Accuracy")+scale_fill_brewer(palette="Dark2")

df2 <- dplyr::filter(df, field == 2)
p2<- ggplot(data=df2, aes(x=depth, y=accuracy, color=dataset)) +
  geom_bar(stat="identity", position=position_dodge())+ #"#E69F00"
  # theme_minimal()+
  theme_classic()+
  geom_line()+
  geom_point()+
  xlab("Depth")+
  ylim(c(0.8,1))+
   ggtitle('2-field')+
  ylab("Accuracy")+scale_fill_brewer(palette="Dark2")


df3 <- dplyr::filter(df, field == 3)
  p3<- ggplot(data=df3, aes(x=depth, y=accuracy, color=dataset)) +
  geom_bar(stat="identity", position=position_dodge())+ #"#E69F00"
  # theme_minimal()+
  theme_classic()+
  geom_line()+
  geom_point()+
  xlab("Depth")+
  ylim(c(0.8,1))+
   ggtitle('3-field')+
  ylab("Accuracy")+scale_fill_brewer(palette="Dark2")

prow <- plot_grid(
  p1 + theme(legend.position="none"),
  p2 + theme(legend.position="none"),
  p3 + theme(legend.position="none"),
  align = 'vh',
  hjust = -1,
  nrow = 1
)
legend <- get_legend(
  p1 +
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom")
)
# legend = get_plot_component(p4, 'guide-box-top', return_all = TRUE)
plot_grid(prow, legend, ncol = 1, rel_heights = c(1, .1))
dev.off()