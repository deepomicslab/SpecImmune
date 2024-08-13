library(ggplot2)

pdf(file="figures/1kg_dp.pdf", width=5, height=3, onefile=FALSE)
df<-read.table("cyp_depth_cutoff.csv", sep=",", header=TRUE)



p <- ggplot(data=df, aes(x=cutoff, y=accuracy, color = method)) +
  geom_line()+
  geom_point()+
  xlab("No. of reads")+
  ylab("Accuracy")+
  scale_color_manual(values = c("slateblue1", "tomato"))

# p + scale_fill_manual(values=c("slateblue1" , "tomato"))
p

dev.off()