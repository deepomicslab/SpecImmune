library(ggplot2)


# pdf(file="figures/hgscv_hifi.pdf", width=5, height=3, onefile=FALSE)
# df<-read.table("HGSCV2_hifi_chain.csv", sep=",", header=TRUE)
# ggplot(data=df, aes(x=cutoff, y=accuracy, color = chain)) +
#   geom_line()+
#   geom_point()+
#   xlab("Depth")+
#   ylab("Accuracy")+theme_classic()+scale_fill_brewer(palette="Dark2")
# dev.off()

pdf(file="figures/HPRC_hifi.pdf", width=5, height=3, onefile=FALSE)
df<-read.table("HPRC_hifi_chain.csv", sep=",", header=TRUE)
ggplot(data=df, aes(x=cutoff, y=accuracy, color = chain)) +
  geom_line()+
  geom_point()+
  xlab("Depth")+
  ylab("Accuracy")+theme_classic()+scale_fill_brewer(palette="Dark2")
dev.off()


