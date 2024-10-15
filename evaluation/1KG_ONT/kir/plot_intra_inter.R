library(ggplot2)
library(ggpubr)

df <- read.table("intra_inter_LD.csv", sep=",", header=TRUE)
pdf(file="intra_inter_LD.pdf", width=4, height=10, onefile=FALSE)
head(df)

vysg <- ggplot(df, aes(x=type,y=Wn,fill=type)) + 
          geom_boxplot() + 
        theme_bw()+
  scale_fill_manual(values = c("#9fc3d5", "#2a347a"))
vysg<-vysg+facet_wrap(Class ~  ., ncol = 1)+
xlab('')+
  theme(legend.position = "none")

vysg+ stat_compare_means(method = "t.test")

dev.off()