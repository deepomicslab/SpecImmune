## load the csv file as a matrix
library(ggplot2)
library(aplot)


# pdf(file="figures/LD_heatmap.pdf", width=6, height=5, onefile=FALSE)
# df<-read.table("kir_LD_values.csv", sep=",", header=TRUE)


# pdf(file="figures/hla_LD_heatmap.pdf", width=6, height=5, onefile=FALSE)
# df<-read.table("hla_LD_values.csv", sep=",", header=TRUE)

pdf(file="figures/hla_kir_LD_heatmap.pdf", width=9, height=8, onefile=FALSE)
df<-read.table("hla_kir_LD_values.csv", sep=",", header=TRUE)

# df <- df[df$sample_num >100,]
## filter NA
df <- df[!is.na(df$Wn),]


# Create a matrix for clustering
df_matrix <- reshape2::acast(df, gene1 ~ gene2, value.var = "Wn", fill = 0)
# df_matrix <- reshape2::acast(df, gene1 ~ gene2, value.var = "D", fill = 0)

# Compute the distance matrix
dist_matrix <- dist(df_matrix)

# Perform hierarchical clustering
hc <- hclust(dist_matrix)

# Reorder the data based on clustering
df$gene1 <- factor(df$gene1, levels = hc$labels[hc$order])
df$gene2 <- factor(df$gene2, levels = hc$labels[hc$order])


# Filter the data to include only the lower triangular part
# df <- df[df$gene1 <= df$gene2, ]

p1<-ggplot(data = df, aes(x=gene1, y=gene2, fill=Wn)) + 
  geom_tile()+
   theme(axis.text.x = element_text())+  #size = 8
   scale_x_discrete(guide = guide_axis(angle = 90))+
 geom_tile(color = "white")+
 scale_fill_gradient2(low = "blue", high = "red", mid = "yellow", midpoint = 0.6,
    space = "Lab",   #midpoint = 0.3, limit = c(0,0.5),
    name="Wn") +
#   theme_classic()

    theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank())

p1

df2<-read.csv("gene_type.csv",header=T, sep=",")
df2$y <- factor(df2$y, levels = hc$labels[hc$order])

p2<-ggplot(df2,aes(x=x,y=y))+
  geom_tile(aes(fill=Type))+
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