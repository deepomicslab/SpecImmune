## load the csv file as a matrix
library(ggplot2)
library(aplot)

library(ggplot2)
library(aplot)
library(reshape2)
library(ggdendro)
library(cowplot)
library(ape)


# pdf(file="figures/hla_kir_LD_heatmap.pdf", width=9, height=8, onefile=FALSE)
# df<-read.table("hla_kir_LD_values.csv", sep=",", header=TRUE)
# df2<-read.csv("gene_type.csv",header=T, sep=",")
# tree_file <- "tree/hla_kir_tree.nwk"

# pdf(file="figures/hla_kir_cyp_vdj_LD_heatmap.pdf", width=25, height=22, onefile=FALSE)
# df<-read.table("hla_kir_cyp_vdj_LD_values.csv", sep=",", header=TRUE)
# df2<-read.csv("vdj_gene_type.csv",header=T, sep=",")
# tree_file <- "tree/all_tree.nwk"

# pdf(file="figures/hla_kir_cyp_vdj_LD_heatmap2.pdf", width=25, height=22, onefile=FALSE)
# df<-read.table("hla_kir_cyp_vdj_LD_values2.csv", sep=",", header=TRUE)
# df2<-read.csv("vdj_gene_type2.csv",header=T, sep=",")
# tree_file <- "tree/all_tree2.nwk"

pdf(file="figures/hla_kir_cyp_LD_heatmap.pdf", width=9, height=8, onefile=FALSE)
df<-read.table("hla_kir_cyp_LD_values.csv", sep=",", header=TRUE)
df2<-read.csv("vdj_gene_type3.csv",header=T, sep=",")
tree_file <- "tree/all_cyp_tree.nwk"

### remove the elements with Wn is NAN or Inf, use filter
df <- df[df$Wn != "Inf",]
df <- df[df$Wn != "NaN",]
# df <- df[df$sample_num >300,]
## filter NA
df <- df[!is.na(df$Wn),]

## remove the df with gene1 start with IG or TCR
df_min <- df[!grepl("^IG", df$gene1),]
df_min <- df_min[!grepl("^TR", df_min$gene1),]
df_min <- df_min[!grepl("^IG", df_min$gene2),]
df_min <- df_min[!grepl("^TR", df_min$gene2),]
## output the df_min to csv file
write.csv(df_min, "hla_kir_cyp_LD_values.csv")


head(df)

## check if IGHD3-10 in df
# df[df$gene1 == "IGHD3-10" | df$gene2 == "IGHD3-10",]

### output the df with Wn is 1
# df1 <- df[df$Wn == 1,]
# head(df1)

# df3 <- df[df$min_w > 0.2,]
# ## covert df3 to a matrix

# df_matrix3 <- reshape2::acast(df3, gene1 ~ gene2, value.var = "min_w", fill = "min_w")
# ### save the matrix to csv file
# write.csv(df_matrix3, "minALD_matrix.tsv")

# Create a matrix for clustering
df_matrix <- reshape2::acast(df, gene1 ~ gene2, value.var = "Wn", fill = 0)

# Compute the distance matrix
dist_matrix <- dist(df_matrix)

# Perform hierarchical clustering
hc <- hclust(dist_matrix, method = "complete")

# Output the tree in nwk format
tree <- as.phylo(hc)
write.tree(tree, file = tree_file)

# Reorder the data based on clustering
df$gene1 <- factor(df$gene1, levels = hc$labels[hc$order])
df$gene2 <- factor(df$gene2, levels = hc$labels[hc$order])


# Filter the data to include only the lower triangular part
# df <- df[df$gene1 <= df$gene2, ]

p1<-ggplot(data = df, aes(x=gene1, y=gene2, fill=Wab)) + 
  geom_tile()+
   theme(axis.text.x = element_text())+  #size = 8
   scale_x_discrete(guide = guide_axis(angle = 90))+
 geom_tile(color = "white")+
scale_fill_gradientn(colors = c("#7eced6", "#d4edee", "#f4dedc", "#f7b8b7", "#f39289"), values = c(0, 0.25, 0.5, 0.75, 1),
# scale_fill_gradientn(colors = c("#7eced6", "#d4edee", "#f4dedc", "#f7b8b7", "#f39289"), values = c(0, 0.1, 0.2, 0.6, 1),
       space = "Lab", name = "Wa/b") +
    theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank())

p1


## remove the df2 with y not in hc$order
df2 <- df2[df2$y %in% hc$labels[hc$order],]
df2$y <- factor(df2$y, levels = hc$labels[hc$order])
# head(hc$labels)
# head(hc$labels[hc$order])
# head(df2$y)

p2<-ggplot(df2,aes(x=x,y=y))+
  geom_tile(aes(fill=Type))+
  scale_x_continuous(expand = c(0,0))+
  theme(panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = "left",
        legend.title = element_blank())+
  scale_fill_manual(values = c("#d9e6eb", "#9fc3d5", "#8f96bd", "#2a347a", "#d6d69b"))





p1%>%
  insert_left(p2,width = 0.05)

# dendro_plot

dev.off()


pdf(file="figures/hla_kir_cyp_LD_heatmap2.pdf", width=9, height=8, onefile=FALSE)
p1<-ggplot(data = df, aes(x=gene1, y=gene2, fill=Wn)) + 
  geom_tile()+
   theme(axis.text.x = element_text())+  #size = 8
   scale_x_discrete(guide = guide_axis(angle = 90))+
 geom_tile(color = "white")+
scale_fill_gradientn(colors = c("#7eced6", "#d4edee", "#f4dedc", "#f7b8b7", "#f39289"), values = c(0, 0.25, 0.5, 0.75, 1),
       space = "Lab", name = "Wn") +
    theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank())

p1


## remove the df2 with y not in hc$order
df2 <- df2[df2$y %in% hc$labels[hc$order],]
df2$y <- factor(df2$y, levels = hc$labels[hc$order])
# head(hc$labels)
# head(hc$labels[hc$order])
# head(df2$y)

p2<-ggplot(df2,aes(x=x,y=y))+
  geom_tile(aes(fill=Type))+
  scale_x_continuous(expand = c(0,0))+
  theme(panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = "left",
        legend.title = element_blank())+
  scale_fill_manual(values = c("#d9e6eb", "#9fc3d5", "#8f96bd", "#2a347a", "#d6d69b"))





p1%>%
  insert_left(p2,width = 0.05)

# dendro_plot

dev.off()


# pdf(file="figures/hla_kir_cyp_vdj_LD_tree.pdf", width=10, height=25, onefile=FALSE)
# # Create the dendrogram plot
# dendro_data <- as.dendrogram(hc)
# dendro_plot <- ggdendrogram(dendro_data, rotate = FALSE) +
#   theme(
#         axis.title.x = element_blank(),
#         axis.title.y = element_blank())

# p1%>%
#   insert_left(dendro_plot,width = 0.05)
# dev.off()