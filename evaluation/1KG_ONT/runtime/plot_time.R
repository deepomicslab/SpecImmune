library(ggplot2)


pdf(file="time_dot.pdf", width=10, height=8, onefile=FALSE)
df<-read.table("vdj_time.csv", sep=",", header=TRUE)

# ggplot(df, aes(x=EUR, y=EAS)) + geom_point()
# Change the point size, and shape
ggplot(df, aes(x=EUR, y=EAS, color=locus)) +
  geom_point(size=1, shape=1)+ 
  geom_text(label=df$Allele,size=2)

dev.off()
