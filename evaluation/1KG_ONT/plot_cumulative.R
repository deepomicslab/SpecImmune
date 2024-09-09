library(ggplot2)


pdf(file="hla/hla_cumulative.pdf", width=6, height=4, onefile=FALSE)
df<-read.table("hla/hla_cumulative.csv", sep=",", header=TRUE)

# Stacked barplot with multiple groups

ggplot(data=df, aes(x=Sample_Num, y=Allele_Num))+
  geom_line()

dev.off()

