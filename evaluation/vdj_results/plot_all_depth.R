library(ggplot2)
library(grid)
library(cowplot)

pdf(file="figures/all_depth.pdf", width=4, height=2, onefile=FALSE)
df<-read.table("all_loci_depth.csv", sep=",", header=TRUE)
ggplot(data=df, aes(x=depth, y=accuracy, color=dataset)) +
  geom_line(size=1)+
  geom_point()+
  xlab("Depth")+
  ylim(c(0.95,1))+
  geom_line(aes(linetype=dataset))+
  ylab("Accuracy")+theme_classic()+ scale_color_manual(values = c("#d9e6eb", "#9fc3d5", "#8f96bd", "#2a347a", "#d6d69b"))+#  scale_color_brewer(palette="Dark2")+
  theme(legend.position="right")
dev.off()
#   geom_line(aes(linetype=dataset))+
#   geom_point(aes(shape=dataset))+
