library(ggplot2)
library(grid)
library(cowplot)

pdf(file="figures/all_depth.pdf", width=6, height=3, onefile=FALSE)
df<-read.table("all_loci_depth.csv", sep=",", header=TRUE)
ggplot(data=df, aes(x=depth, y=accuracy, color=dataset)) +
  geom_line()+
  geom_point()+
  xlab("Depth")+
  ylim(c(0.95,1))+
  ylab("Accuracy")+theme_classic()+scale_color_brewer(palette="Dark2")+
  theme(legend.position="none")
dev.off()
#   geom_line(aes(linetype=dataset))+
#   geom_point(aes(shape=dataset))+
