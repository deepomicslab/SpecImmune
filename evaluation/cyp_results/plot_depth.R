library(ggplot2)


# pdf(file="figures/1kg_dp.pdf", width=5, height=3, onefile=FALSE)
# df<-read.table("cyp_depth_cutoff.csv", sep=",", header=TRUE)
# ggplot(data=df, aes(x=cutoff, y=accuracy, color = method)) +
#   geom_line()+
#   geom_point()+
#   xlab("No. of reads")+
#   ylab("Accuracy")+
#   scale_color_manual(values = c("#9fc3d5","#2a347a"))+theme_classic()
# dev.off()

pdf(file="figures/hprc_hifi.pdf", width=5, height=3, onefile=FALSE)
df<-read.table("hprc_hifi_cyp_depth_cutoff.csv", sep=",", header=TRUE)
ggplot(data=df, aes(x=cutoff, y=accuracy, color = method)) +
  geom_line()+
  geom_point()+
  xlab("No. of reads")+
  ylab("Accuracy")+
  scale_color_manual(values = c("#9fc3d5","#2a347a"))+theme_classic()
dev.off()

pdf(file="figures/hprc_ont.pdf", width=5, height=3, onefile=FALSE)
df<-read.table("hprc_ont_cyp_depth_cutoff.csv", sep=",", header=TRUE)
ggplot(data=df, aes(x=cutoff, y=accuracy, color = method)) +
  geom_line()+
  geom_point()+
  xlab("No. of reads")+
  ylab("Accuracy")+
  scale_color_manual(values = c("#9fc3d5","#2a347a"))+theme_classic()
dev.off()

