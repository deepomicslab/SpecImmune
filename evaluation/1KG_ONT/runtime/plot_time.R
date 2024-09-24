library(ggplot2)


pdf(file="time_dot.pdf", width=5, height=4, onefile=FALSE)
df<-read.table("vdj_time.csv", sep=",", header=TRUE)

# ggplot(df, aes(x=EUR, y=EAS)) + geom_point()
# Change the point size, and shape
ggplot(df, aes(x=time, y=mem)) +
  geom_point()+
  xlab("CPU Time (h)") +
  ylab("Peak Memory (GB)")+
  ggtitle("IG+TCR")
  # + 
  # geom_text(label=df$Allele,size=2)
## count mean value of time and mem
print(paste("mean time: ", mean(df$time)))
print(paste("mean mem: ", mean(df$mem)))

dev.off()
