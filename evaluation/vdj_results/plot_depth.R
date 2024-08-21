library(ggplot2)


# pdf(file="figures/1kg_dp.pdf", width=5, height=3, onefile=FALSE)
# df<-read.table("cyp_depth_cutoff.csv", sep=",", header=TRUE)
# ggplot(data=df, aes(x=cutoff, y=accuracy, color = method)) +
#   geom_line()+
#   geom_point()+
#   xlab("No. of reads")+
#   ylab("Accuracy")+
#   scale_color_manual(values = c("slateblue1", "tomato"))+theme_classic()
# dev.off()

pdf(file="figures/hgscv_hifi.pdf", width=5, height=3, onefile=FALSE)
df<-read.table("HGSCV2_hifi_chain.csv", sep=",", header=TRUE)
ggplot(data=df, aes(x=cutoff, y=accuracy, color = chain)) +
  geom_line()+
  geom_point()+
  xlab("No. of reads")+
  ylab("Accuracy")+theme_classic()+scale_fill_brewer(palette="Dark2")
dev.off()


