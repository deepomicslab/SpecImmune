library(ggplot2)


pdf(file="hla/hla_dot.pdf", width=10, height=8, onefile=FALSE)
df<-read.table("hla/hla_freq.csv", sep=",", header=TRUE)

# ggplot(df, aes(x=EUR, y=EAS)) + geom_point()
# Change the point size, and shape
ggplot(df, aes(x=EUR, y=EAS, color=locus)) +
  geom_point(size=1, shape=1)

dev.off()