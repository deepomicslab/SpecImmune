library(ggplot2)
library(grid)
library(cowplot)
library(dplyr)


pdf(file="figures/kir_depth_change2.pdf", width=8, height=4, onefile=FALSE)

df<-read.table("all_loci_depth.csv", sep=",", header=TRUE)
# df <- dplyr::filter(df, depth != 20)
df <- dplyr::filter(df, depth == 0)
df
## make accueacy as numeric
df$accuracy <- as.numeric(as.character(df$accuracy))

# p1<- ggplot(data=df, aes(x=dataset, y=accuracy)) +
#   geom_bar(stat="identity", position=position_dodge(), aes(fill=dataset))+ #"#E69F00"
#   ylab("Accuracy")+ylim(c(0.5,1))+scale_fill_brewer(palette="Dark2")+ 
#   geom_tile()+
#    theme(axis.text.x = element_text(size = 8))+
#    scale_x_discrete(guide = guide_axis(angle = 90))
df1 <- dplyr::filter(df, field == 1)

p1<-ggplot(df1, aes(x=dataset, y=accuracy, fill=dataset)) +
  geom_bar(stat="identity") +
  theme_classic() +
  scale_fill_manual(values = c("#d9e6eb", "#9fc3d5", "#8f96bd", "#2a347a", "#d6d69b")) +
  coord_cartesian(ylim = c(0.85, 0.95)) +
  theme(axis.text.x = element_text(size = 8)) +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  ggtitle('1-field') +
  xlab('')

df2 <- dplyr::filter(df, field == 2)

p2<-ggplot(df2, aes(x=dataset, y=accuracy, fill=dataset)) +
  geom_bar(stat="identity")+theme_classic()+scale_fill_manual(values = c("#d9e6eb", "#9fc3d5", "#8f96bd", "#2a347a", "#d6d69b"))+ coord_cartesian(ylim=c(0.85,0.95))+ 
   theme(axis.text.x = element_text(size = 8))+
   scale_x_discrete(guide = guide_axis(angle = 90))+ggtitle('2-field')+xlab('')


df3 <- dplyr::filter(df, field == 3)

p3<-ggplot(df3, aes(x=dataset, y=accuracy, fill=dataset)) +
  geom_bar(stat="identity")+theme_classic()+scale_fill_manual(values = c("#d9e6eb", "#9fc3d5", "#8f96bd", "#2a347a", "#d6d69b"))+ coord_cartesian(ylim=c(0.85,0.95))+ 
   theme(axis.text.x = element_text(size = 8))+
   scale_x_discrete(guide = guide_axis(angle = 90))+ggtitle('3-field')+xlab('')


prow <- plot_grid(
  p1 + theme(legend.position="none"),
  p2 + theme(legend.position="none"),
  p3 + theme(legend.position="none"),
  align = 'vh',
  hjust = -1,
  nrow = 1
)
legend <- get_legend(
  p1 +
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom")
)
df1
# legend = get_plot_component(p4, 'guide-box-top', return_all = TRUE)
plot_grid(prow, legend, ncol = 1, rel_heights = c(1, .1))

dev.off()

