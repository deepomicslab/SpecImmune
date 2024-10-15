library(ggplot2)
library(grid)
library(cowplot)
library(dplyr)


pdf(file="figures/cyp_depth.pdf", width=10, height=3, onefile=FALSE)


df<-read.table("cyp_depth_cutoff.csv", sep=",", header=TRUE)
p1<-ggplot(data=df, aes(x=cutoff, y=accuracy, color = method)) +
  geom_line()+
  geom_point()+
  xlab("No. of reads")+
  ylab("Accuracy")+
  ggtitle('1kGP') +
  scale_color_manual(values = c("#9fc3d5","#2a347a"))+theme_classic()



df<-read.table("hprc_hifi_cyp_depth_cutoff.csv", sep=",", header=TRUE)
p2<-ggplot(data=df, aes(x=cutoff, y=accuracy, color = method)) +
  geom_line()+
  geom_point()+
  xlab("No. of reads")+
  ylab("Accuracy")+
  ggtitle('HPRC HiFi') +
  scale_color_manual(values = c("#9fc3d5","#2a347a"))+theme_classic()



df<-read.table("hprc_ont_cyp_depth_cutoff.csv", sep=",", header=TRUE)
p3<-ggplot(data=df, aes(x=cutoff, y=accuracy, color = method)) +
  geom_line()+
  geom_point()+
  xlab("No. of reads")+
  ylab("Accuracy")+
  ggtitle('HPRC ONT') +
  scale_color_manual(values = c("#9fc3d5","#2a347a"))+theme_classic()


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

