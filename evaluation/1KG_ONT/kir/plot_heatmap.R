## load the csv file as a matrix
library(ggplot2)
library(aplot)


pdf(file="figures/fst_heatmap.pdf", width=6, height=5, onefile=FALSE)


df<-read.table("kir_fst.csv", sep=",", header=TRUE)

## assign x order
df$Group1 <- factor(df$Group1, levels = c('CDX', 'CHB', 'CHS', 'JPT', 'KHV', 'BEB', 'GIH', 'ITU', 'PJL', 'STU', 'ACB', 'ASW', 'ESN', 'GWD', 'LWK', 'MSL', 'YRI', 'CEU', 'FIN', 'GBR', 'IBS', 'TSI', 'CLM', 'MXL', 'PEL', 'PUR'))
df$Group2 <- factor(df$Group2, levels = c('CDX', 'CHB', 'CHS', 'JPT', 'KHV', 'BEB', 'GIH', 'ITU', 'PJL', 'STU', 'ACB', 'ASW', 'ESN', 'GWD', 'LWK', 'MSL', 'YRI', 'CEU', 'FIN', 'GBR', 'IBS', 'TSI', 'CLM', 'MXL', 'PEL', 'PUR'))


# Create a matrix for clustering
df_matrix <- reshape2::acast(df, Group1 ~ Group2, value.var = "JSD", fill = 0)

# Compute the distance matrix
dist_matrix <- dist(df_matrix)

# Perform hierarchical clustering
hc <- hclust(dist_matrix)

# Reorder the data based on clustering
df$Group1 <- factor(df$Group1, levels = hc$labels[hc$order])
df$Group2 <- factor(df$Group2, levels = hc$labels[hc$order])


p1<-ggplot(data = df, aes(x=Group1, y=Group2, fill=JSD)) + 
  geom_tile()+
   theme(axis.text.x = element_text(size = 8))+
   scale_x_discrete(guide = guide_axis(angle = 90))+
 geom_tile(color = "white")+
#  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0.6,
scale_fill_gradient2(low = "#8fb4be", high = "#d93f49", mid = "white", midpoint = 0.5,
# scale_fill_gradientn(colors = c("#d93f49", "#e28187", "#ebbfc2", "#d5e1e3", "#afc9cf","#8fb4be"), values = c(0, 0.3, 0.6, 0.62,0.64, 0.66),
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

df2<-read.csv("kir_color.csv",header=T, sep=",")
# df2$y<-factor(df2$y,levels = rev(df2$y))
df2$y <- factor(df2$y, levels = hc$labels[hc$order])


p2<-ggplot(df2,aes(x=x,y=y))+
  geom_tile(aes(fill=group))+
  scale_x_continuous(expand = c(0,0))+
  theme(panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = "left",
        legend.title = element_blank())+
  scale_fill_manual(values = c("#d9e6eb", "#9fc3d5", "#8f96bd", "#2a347a", "#d6d69b"))
#         +
#   scale_fill_manual(values = c("green","blue","red"))
p1%>%
  insert_left(p2,width = 0.05)

## remove the df with compare2 is "same"
df3 <- df[df$Group1 != df$Group2,]
## compare the value with compare of intra and inter using t-test
ttest <- t.test(df3$JSD[df3$compare == "intra"], df3$JSD[df3$compare == "inter"])
ttest


dev.off()