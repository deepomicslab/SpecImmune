library(ggplot2)
library(grid)
library(cowplot)
library(ggrepel)


pdf(file="figures/vdj_depth.pdf", width=7, height=6, onefile=FALSE)

df<-read.table("HGSCV2_hifi_chain.csv", sep=",", header=TRUE)
# Only label endpoints (depth 0 and 20)
df$label <- ifelse(df$cutoff %in% c(0, 20), df$total, "")
p1<- ggplot(data=df, aes(x=cutoff, y=accuracy, group = chain)) +
  geom_line(aes(color=chain), size=1)+
  geom_point(aes(color=chain))+
  geom_text_repel(aes(label=label, color=chain), size=2.5, show.legend=FALSE, 
                  max.overlaps = 50, box.padding = 0.5, point.padding = 0.3,
                  min.segment.length = 0.1, force = 3, seed = 42, 
                  segment.size = 0.3, segment.alpha = 0.5) +
  xlab("Minimum depth (filter cutoff)")+
  ylim(c(0.85,1.02))+
  ggtitle('HGSVC HiFi')+
  ylab("Accuracy")+theme_classic()+
  scale_color_manual(values = c("#827e3f", "skyblue", "#d9e6eb", "#9fc3d5", "#8f96bd", "#2a347a", "#d6d69b")) 


df<-read.table("hgscv2_clr_chain.csv", sep=",", header=TRUE)
df$label <- ifelse(df$cutoff %in% c(0, 20), df$total, "")
p2<-ggplot(data=df, aes(x=cutoff, y=accuracy, color = chain)) +
  geom_line(size=1)+
  geom_point()+
  geom_text_repel(aes(label=label, color=chain), size=2.5, show.legend=FALSE, 
                  max.overlaps = 50, box.padding = 0.5, point.padding = 0.3,
                  min.segment.length = 0.1, force = 3, seed = 42, 
                  segment.size = 0.3, segment.alpha = 0.5) +
  xlab("Minimum depth (filter cutoff)")+
  ylim(c(0.85,1.02))+
  ggtitle('HGSVC CLR')+
  ylab("Accuracy")+theme_classic()+
  scale_color_manual(values = c("#827e3f", "skyblue", "#d9e6eb", "#9fc3d5", "#8f96bd", "#2a347a", "#d6d69b")) 

df<-read.table("HPRC_hifi_chain.csv", sep=",", header=TRUE)
df$label <- ifelse(df$cutoff %in% c(0, 20), df$total, "")
p3<-ggplot(data=df, aes(x=cutoff, y=accuracy, color = chain)) +
  geom_line(size=1)+
  geom_point()+
  geom_text_repel(aes(label=label, color=chain), size=2.5, show.legend=FALSE, 
                  max.overlaps = 50, box.padding = 0.5, point.padding = 0.3,
                  min.segment.length = 0.1, force = 3, seed = 42, 
                  segment.size = 0.3, segment.alpha = 0.5) +
  xlab("Minimum depth (filter cutoff)")+
  ylim(c(0.85,1.02))+
  ggtitle('HPRC HiFi')+
  ylab("Accuracy")+theme_classic()+
  scale_color_manual(values = c("#827e3f", "skyblue", "#d9e6eb", "#9fc3d5", "#8f96bd", "#2a347a", "#d6d69b")) 

df<-read.table("HPRC_ont_chain.csv", sep=",", header=TRUE)
df$label <- ifelse(df$cutoff %in% c(0, 20), df$total, "")
p4<-ggplot(data=df, aes(x=cutoff, y=accuracy, color = chain)) +
  geom_line(size=1)+
  geom_point()+
  geom_text_repel(aes(label=label, color=chain), size=2.5, show.legend=FALSE, 
                  max.overlaps = 50, box.padding = 0.5, point.padding = 0.3,
                  min.segment.length = 0.1, force = 3, seed = 42, 
                  segment.size = 0.3, segment.alpha = 0.5) +
  xlab("Minimum depth (filter cutoff)")+
  ylim(c(0.85,1.02))+
  ggtitle('HPRC ONT')+
  ylab("Accuracy")+theme_classic()+
  scale_color_manual(values = c("#827e3f", "skyblue", "#d9e6eb", "#9fc3d5", "#8f96bd", "#2a347a", "#d6d69b")) 

prow <- plot_grid(
  p3 + theme(legend.position="none"),
  p4 + theme(legend.position="none"),
  p1 + theme(legend.position="none"),
  p2 + theme(legend.position="none"),
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

