library(ggplot2)
library(grid)
library(cowplot)


pdf(file="figures/vdj_depth.pdf", width=6, height=5, onefile=FALSE)

df<-read.table("HGSCV2_hifi_chain.csv", sep=",", header=TRUE)
p1<- ggplot(data=df, aes(x=cutoff, y=accuracy, group = chain)) +
  geom_line(aes(color=chain), size=1)+
  geom_point(aes(color=chain))+
  xlab("Depth")+
  ylim(c(0.85,1))+
  ggtitle('HGSVC HiFi')+
  ylab("Accuracy")+theme_classic()+
  scale_color_manual(values = c("#827e3f", "skyblue", "#d9e6eb", "#9fc3d5", "#8f96bd", "#2a347a", "#d6d69b")) 
  # +scale_color_brewer(palette="Dark2")


df<-read.table("hgscv2_clr_chain.csv", sep=",", header=TRUE)
p2<-ggplot(data=df, aes(x=cutoff, y=accuracy, color = chain)) +
  geom_line(size=1)+
  geom_point()+
  xlab("Depth")+
  ylim(c(0.85,1))+
  ggtitle('HGSVC CLR')+
  ylab("Accuracy")+theme_classic()+
  scale_color_manual(values = c("#827e3f", "skyblue", "#d9e6eb", "#9fc3d5", "#8f96bd", "#2a347a", "#d6d69b")) 

df<-read.table("HPRC_hifi_chain.csv", sep=",", header=TRUE)
p3<-ggplot(data=df, aes(x=cutoff, y=accuracy, color = chain)) +
  geom_line(size=1)+
  geom_point()+
  xlab("Depth")+
  ylim(c(0.85,1))+
  ggtitle('HPRC HiFi')+
  ylab("Accuracy")+theme_classic()+
  scale_color_manual(values = c("#827e3f", "skyblue", "#d9e6eb", "#9fc3d5", "#8f96bd", "#2a347a", "#d6d69b")) 

df<-read.table("HPRC_ont_chain.csv", sep=",", header=TRUE)
p4<-ggplot(data=df, aes(x=cutoff, y=accuracy, color = chain)) +
  geom_line(size=1)+
  geom_point()+
  xlab("Depth")+
  ylim(c(0.85,1))+
  ggtitle('HPRC ONT')+
  ylab("Accuracy")+theme_classic()+
  scale_color_manual(values = c("#827e3f", "skyblue", "#d9e6eb", "#9fc3d5", "#8f96bd", "#2a347a", "#d6d69b")) 

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
  p4 +
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom")
)
# legend = get_plot_component(p4, 'guide-box-top', return_all = TRUE)
plot_grid(prow, legend, ncol = 1, rel_heights = c(1, .1))



dev.off()



pdf(file="figures/hgscv_hifi.pdf", width=5, height=3, onefile=FALSE)
df<-read.table("HGSCV2_hifi_chain.csv", sep=",", header=TRUE)
ggplot(data=df, aes(x=cutoff, y=accuracy, color = chain)) +
  geom_line()+
  geom_point()+
  xlab("Depth")+
  ylab("Accuracy")+theme_classic()+
  scale_color_manual(values = c("#827e3f", "skyblue", "#d9e6eb", "#9fc3d5", "#8f96bd", "#2a347a", "#d6d69b"))
dev.off()

pdf(file="figures/HGSCV_CLR.pdf", width=5, height=3, onefile=FALSE)
df<-read.table("hgscv2_clr_chain.csv", sep=",", header=TRUE)
ggplot(data=df, aes(x=cutoff, y=accuracy, color = chain)) +
  geom_line()+
  geom_point()+
  xlab("Depth")+
  ylab("Accuracy")+theme_classic()+
  scale_color_manual(values = c("#827e3f", "skyblue", "#d9e6eb", "#9fc3d5", "#8f96bd", "#2a347a", "#d6d69b"))
dev.off()

pdf(file="figures/HPRC_ont.pdf", width=5, height=3, onefile=FALSE)
df<-read.table("HPRC_ont_chain.csv", sep=",", header=TRUE)
ggplot(data=df, aes(x=cutoff, y=accuracy, color = chain)) +
  geom_line()+
  geom_point()+
  xlab("Depth")+
  ylab("Accuracy")+theme_classic()+
  scale_color_manual(values = c("#827e3f", "skyblue", "#d9e6eb", "#9fc3d5", "#8f96bd", "#2a347a", "#d6d69b"))
dev.off()

pdf(file="figures/HPRC_hifi.pdf", width=5, height=3, onefile=FALSE)
df<-read.table("HPRC_hifi_chain.csv", sep=",", header=TRUE)
ggplot(data=df, aes(x=cutoff, y=accuracy, color = chain)) +
  geom_line()+
  geom_point()+
  xlab("Depth")+
  ylab("Accuracy")+theme_classic()+
  scale_color_manual(values = c("#827e3f", "skyblue","#d9e6eb", "#9fc3d5", "#8f96bd", "#2a347a", "#d6d69b"))
dev.off()

