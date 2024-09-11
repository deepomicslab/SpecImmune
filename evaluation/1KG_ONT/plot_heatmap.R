## load the csv file as a matrix
library(ggplot2)
library(aplot)


pdf(file="hla/fst_heatmap.pdf", width=6, height=5, onefile=FALSE)


df<-read.table("hla/hla_fst.csv", sep=",", header=TRUE)

## assign x order
df$Group1 <- factor(df$Group1, levels = c('CDX', 'CHB', 'CHS', 'JPT', 'KHV', 'BEB', 'GIH', 'ITU', 'PJL', 'STU', 'ACB', 'ASW', 'ESN', 'GWD', 'LWK', 'MSL', 'YRI', 'CEU', 'FIN', 'GBR', 'IBS', 'TSI', 'CLM', 'MXL', 'PEL', 'PUR'))
df$Group2 <- factor(df$Group2, levels = c('CDX', 'CHB', 'CHS', 'JPT', 'KHV', 'BEB', 'GIH', 'ITU', 'PJL', 'STU', 'ACB', 'ASW', 'ESN', 'GWD', 'LWK', 'MSL', 'YRI', 'CEU', 'FIN', 'GBR', 'IBS', 'TSI', 'CLM', 'MXL', 'PEL', 'PUR'))



p1<-ggplot(data = df, aes(x=Group1, y=Group2, fill=JSD)) + 
  geom_tile()+
   theme(axis.text.x = element_text(size = 8))+
   scale_x_discrete(guide = guide_axis(angle = 90))+
 geom_tile(color = "white")+
 scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0.4,
    space = "Lab",   #midpoint = 0.3, limit = c(0,0.5),
    name="JSD") +
#   theme_classic()

    theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank())

df2<-read.csv("hla/hla_color.csv",header=T, sep=",")
df2$y<-factor(df2$y,levels = rev(df2$y))
p2<-ggplot(df2,aes(x=x,y=y))+
  geom_tile(aes(fill=group))+
  scale_x_continuous(expand = c(0,0))+
  theme(panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = "left",
        legend.title = element_blank())
#         +
#   scale_fill_manual(values = c("green","blue","red"))
p1%>%
  insert_left(p2,width = 0.05)

dev.off()