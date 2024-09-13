## load the csv file as a matrix
library(ggplot2)
library(cowplot)
library(dplyr)


df<-read.table("cyp_activity.csv", sep=",", header=TRUE)

  pdf(file="cyp_activity_box.pdf", width=3, height=3, onefile=FALSE)
  
  p1<-ggplot(df, aes(x=Super_Pop, y = Activity_score, fill=Super_Pop)) +
    geom_boxplot()  +
    stat_summary(fun=mean, geom="point", shape=20, size=5, color="yellow", fill="yellow") +
    theme_classic()+
    xlab('')+
    ylab('CYP2D6 Activity Score')+
  theme(legend.position = "none")

  p1+scale_fill_brewer(palette="Dark2")
  

df1<-filter(df, Super_Pop == "EAS")
df2<-filter(df, Super_Pop != "EAS")
## compare the values of the activity score of the two groups using t test

t.test(df1$Activity_score, df2$Activity_score)


  dev.off()