# 读取 PCA 结果
pca_data <- read.table("pca_out.eigenvec", header=FALSE)
colnames(pca_data) <- c("FID", "IID", paste0("PC", 1:(ncol(pca_data)-2)))

# 读取群体信息
pop_info <- read.table("pop.txt", header=FALSE)
colnames(pop_info) <- c("FID", "IID", "Population")

# 合并数据
merged_data <- merge(pca_data, pop_info, by=c("FID", "IID"))

# 绘制 PCA 散点图（以 PC1 和 PC2 为例）
library(ggplot2)
g <- ggplot(merged_data, aes(x=PC1, y=PC2, color=Population)) +
  geom_point(size=2, alpha=0.8) +
  theme_minimal() +
  labs(title="PCA 分析", x="主成分 1 (PC1)", y="主成分 2 (PC2)")

print(g)

library(plotly)

# 绘制三维散点图
g <- plot_ly(merged_data, x = ~PC1, y = ~PC2, z = ~PC3, color = ~Population, colors = c("red", "green", "blue")) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'PC1'),
                      yaxis = list(title = 'PC2'),
                      zaxis = list(title = 'PC3')),
         title = 'PCA分析 - 三维散点图')

print(g)