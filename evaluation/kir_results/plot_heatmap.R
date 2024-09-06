## load the csv file as a matrix
library(ggplot2)


pdf(file="figures/kir_heatmap.pdf", width=6, height=3, onefile=FALSE)


## row name is true
# df<-read.table("sim_results.csv", sep=",", header=TRUE, row.names=1)
# mat<-as.matrix(df)
# heatmap(mat, Rowv=NA, Colv=NA, col = colorRampPalette(c("white", "blue"))(100), scale="none", margins=c(5,10))


df<-read.table("sim_results2.csv", sep=",", header=TRUE)

ggplot(data = df, aes(x=Gene, y=Depth, fill=Accuracy)) + 
  geom_tile()+
   theme(axis.text.x = element_text(size = 8))+
   scale_x_discrete(guide = guide_axis(angle = 90))+
 geom_tile(color = "white")+
 scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
   midpoint = 0.8, limit = c(0.65,1), space = "Lab", 
    name="3-field\nAccuracy") +
#   theme_classic()

    theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank())


dev.off()