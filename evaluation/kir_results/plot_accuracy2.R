library(ggplot2)
library(grid)
library(cowplot)
library(dplyr)


pdf(file="figures/kir_depth_change.pdf", width=8, height=3, onefile=FALSE)

df<-read.table("all_loci_depth.csv", sep=",", header=TRUE)
df <- dplyr::filter(df, depth != 20)

## covert the field from numeric to factor
df$field <- as.factor(df$field)

df1 <- dplyr::filter(df, dataset == "HPRC HiFi")
p1<- ggplot(data=df1, aes(x=depth, y=accuracy, group=field)) +
  theme_classic()+
  geom_line(aes(color=field))+
  geom_point(aes(color=field))+
  xlab("No. of reads")+
  # ylim(c(0.85,1))+
  ggtitle('HPRC HiFi')+
  ylab("Accuracy")+scale_color_brewer(palette="Dark2")


df2 <- dplyr::filter(df, dataset == "HPRC ONT")
p2<- ggplot(data=df2, aes(x=depth, y=accuracy, group=field)) +
  theme_classic()+
  geom_line(aes(color=field))+
  geom_point(aes(color=field))+
  xlab("No. of reads")+
  # ylim(c(0.85,1))+
  ggtitle('HPRC ONT')+
  ylab("Accuracy")+scale_color_brewer(palette="Dark2")

df3 <- dplyr::filter(df, dataset == "HGSCV2 HiFi")
p3<- ggplot(data=df3, aes(x=depth, y=accuracy, group=field)) +
  theme_classic()+
  geom_line(aes(color=field))+
  geom_point(aes(color=field))+
  xlab("No. of reads")+
  # ylim(c(0.85,1))+
  ggtitle('HGSCV HiFi')+
  ylab("Accuracy")+scale_color_brewer(palette="Dark2")

df4 <- dplyr::filter(df, dataset == "HGSCV2 CLR")
p4<- ggplot(data=df4, aes(x=depth, y=accuracy, group=field)) +
  theme_classic()+
  geom_line(aes(color=field))+
  geom_point(aes(color=field))+
  xlab("No. of reads")+
  # ylim(c(0.85,1))+
  ggtitle('HGSCV2 CLR')+
  ylab("Accuracy")+scale_color_brewer(palette="Dark2")





prow <- plot_grid(
  p1 + theme(legend.position="none"),
  p2 + theme(legend.position="none"),
  p3 + theme(legend.position="none"),
  p4 + theme(legend.position="none"),
  align = 'vh',
  hjust = -1,
  nrow = 2
)
legend <- get_legend(
  p1 +
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom")
)
df4
# legend = get_plot_component(p4, 'guide-box-top', return_all = TRUE)
plot_grid(prow, legend, ncol = 1, rel_heights = c(1, .1))
dev.off()