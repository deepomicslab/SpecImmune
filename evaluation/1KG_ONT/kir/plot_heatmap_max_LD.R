## load the csv file as a matrix
library(ggplot2)
library(aplot)
library(reshape2)
library(dplyr)


pdf(file="figures/LD_max.pdf", width=3.5, height=2, onefile=FALSE)


df<-read.table("tree/vdj_cyp_hla_kir_max.csv", sep=",", header=TRUE)


  get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
  }
  ## transform df to matrix
  # Transform df to matrix
df_matrix <- acast(df, Class1 ~ Class2, value.var = "ALD")
df_matrix

df <- get_lower_tri(df_matrix)
df
df <- melt(df, na.rm = TRUE)

df


# ggheatmap <- ggplot(df, aes(Class1, Class2, fill = ALD))+
#  geom_tile(color = "white")+
#  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
#    midpoint = 0, limit = c(-1,1), space = "Lab", 
#     name="Pearson\nCorrelation") +
#   theme_minimal()+ # minimal theme
#  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
#     size = 12, hjust = 1))+
#  coord_fixed()
# # Print the heatmap
# print(ggheatmap)

p1<-ggplot(data = df,aes(Var2, Var1, fill = value)) + 
  geom_text(aes(label = "x")) +
  # geom_text(aes(label = round(ALD, 2)), size = 3) +
  #  theme(axis.text.x = element_text(size = 8))+
  #  scale_x_discrete(guide = guide_axis(angle = 90))+
   geom_tile(color = "white")+
#  scale_fill_gradient2(low = "blue", high = "red", mid = "yellow", midpoint = 0.3,
scale_fill_gradient2(low = "#d93f49", high = "#8fb4be",
# scale_fill_gradientn(colors = c("#d93f49", "#e28187", "#ebbfc2", "#d5e1e3", "#afc9cf","#8fb4be"), values = c(0, 0.2, 0.4, 0.6,0.8, 1),
    space = "Lab",   #midpoint = 0.3, limit = c(0,0.5),
    name="Max minALD") +
  # theme_minimal()
    theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank())

p1+geom_text(aes(label = round(value, 2)), size = 4, color="white")


dev.off()