## load the csv file as a matrix
library(ggplot2)
library(aplot)


pdf(file="figures/LD_heatmap.pdf", width=6, height=5, onefile=FALSE)


df<-read.table("kir_LD_values.csv", sep=",", header=TRUE)


p1<-ggplot(data = df, aes(x=gene1, y=gene2, fill=D)) + 
  geom_tile()+
   theme(axis.text.x = element_text(size = 8))+
   scale_x_discrete(guide = guide_axis(angle = 90))+
 geom_tile(color = "white")+
 scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0.6,
    space = "Lab",   #midpoint = 0.3, limit = c(0,0.5),
    name="D") +
#   theme_classic()

    theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank())

p1

dev.off()