library(ggplot2)


pdf(file="hla/hla_cumulative.pdf", width=4, height=3, onefile=FALSE)
df<-read.table("hla/hla_cumulative.csv", sep=",", header=TRUE)

# Stacked barplot with multiple groups

ggplot(data=df, aes(x=Sample_Num, y=Allele_Num))+
  geom_line(color="skyblue",linewidth=1.2)+
  theme_bw()+
  ylab("No. of unique alleles")+xlab("No. of samples")

dev.off()

