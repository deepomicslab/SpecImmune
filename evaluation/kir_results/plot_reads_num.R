library(ggplot2)
library(grid)
library(cowplot)
library(dplyr)

# [1] "#1B9E77" "#D95F02" "#7570B3" "#E7298A" "#66A61E" "#E6AB02" "#A6761D"
# [8] "#666666"

pdf(file="figures/read_num.pdf", width=16, height=6, onefile=FALSE)

df<-read.table("HPRC_ont_read_num.csv", sep=",", header=TRUE)
p1 <- ggplot(data = df, aes(x = gene, y = read_num)) +
  geom_boxplot(fill = "#1B9E77") +
  geom_point(stat = "summary", fun = "mean", shape = 23, size = 3, fill = "white", color = "black") +
  theme_classic() +
  theme(legend.position = "none", axis.text.x = element_text(size = 10)) +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  xlab("") +
  ggtitle('HPRC ONT') +
  ylab("No. of reads") 

df<-read.table("hgscv2_clr_read_num.csv", sep=",", header=TRUE)
p2 <- ggplot(data = df, aes(x = gene, y = read_num)) +
  geom_boxplot(fill = "#1B9E77") +
  geom_point(stat = "summary", fun = "mean", shape = 23, size = 3, fill = "white", color = "black") +
  theme_classic() +
  theme(legend.position = "none", axis.text.x = element_text(size = 10)) +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  xlab("") +
  ggtitle('HGSCV CLR') +
  ylab("No. of reads") 

  df<-read.table("HPRC_hifi_read_num.csv", sep=",", header=TRUE)
p3 <- ggplot(data = df, aes(x = gene, y = read_num)) +
  geom_boxplot(fill = "#1B9E77") +
  geom_point(stat = "summary", fun = "mean", shape = 23, size = 3, fill = "white", color = "black") +
  theme_classic() +
  theme(legend.position = "none", axis.text.x = element_text(size = 10)) +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  xlab("") +
  ggtitle('HPRC HiFi') +
  ylab("No. of reads") 


  df<-read.table("HGSCV2_hifi_read_num.csv", sep=",", header=TRUE)
p4 <- ggplot(data = df, aes(x = gene, y = read_num)) +
  geom_boxplot(fill = "#1B9E77") +
  geom_point(stat = "summary", fun = "mean", shape = 23, size = 3, fill = "white", color = "black") +
  theme_classic() +
  theme(legend.position = "none", axis.text.x = element_text(size = 10)) +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  xlab("") +
  ggtitle('HGSCV HiFi') +
  ylab("No. of reads") 



prow <- plot_grid(
  p1 + theme(legend.position="none"),
  p2 + theme(legend.position="none"),
  p3 + theme(legend.position="none"),
  p4 + theme(legend.position="none"),
  align = 'vh',
  hjust = -1,
  nrow = 2
)
legend <- get_legend(
  p1 +
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom")
)
# legend = get_plot_component(p4, 'guide-box-top', return_all = TRUE)
plot_grid(prow, legend, ncol = 1, rel_heights = c(1, .1))


dev.off()