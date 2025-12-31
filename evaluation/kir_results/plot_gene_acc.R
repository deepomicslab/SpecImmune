library(ggplot2)
library(grid)
library(cowplot)
library(dplyr)


pdf(file="figures/kir_depth.pdf", width=16, height=3, onefile=FALSE)

df<-read.table("HPRC_hifi_field1_result.csv", sep=",", header=TRUE)
# df <- dplyr::filter(df, cutoff == 0)

df2<-read.table("HGSCV2_hifi_field1_result.csv", sep=",", header=TRUE)
# df2 <- dplyr::filter(df2, cutoff == 0)

df3<-read.table("HPRC_ont_field1_result.csv", sep=",", header=TRUE)
# df3 <- dplyr::filter(df3, cutoff == 0)

df4<-read.table("hgscv2_clr_field1_result.csv", sep=",", header=TRUE)
# df4 <- dplyr::filter(df4, cutoff == 0)


# ## add a column to the data frame
# df$dataset <- "HPRC"
# df2$dataset <- "HGSCV2"
# ## combine the two data frames
# df <- rbind(df, df2)
## for all the dfs, only keep rows where cutoff == 0
df <- dplyr::filter(df, cutoff == 0)
df2 <- dplyr::filter(df2, cutoff == 0)
df3 <- dplyr::filter(df3, cutoff == 0)
df4 <- dplyr::filter(df4, cutoff == 0)

# Plot accuracy with count labels on top
p1<- ggplot(data=df, aes(x=gene, y=accuracy)) +
  geom_bar(stat="identity", position=position_dodge(), fill="#2a347a")+ 
  geom_text(aes(label=total), vjust=0.5, hjust=1.5, size=3.5, angle=90, color="white")+
  theme_classic()+
   theme(legend.position = "none", axis.text.x = element_text(size = 10))+
   scale_x_discrete(guide = guide_axis(angle = 90))+
  xlab("")+
  ggtitle('HPRC HiFi')+
  ylab("Accuracy")

p2<- ggplot(data=df2, aes(x=gene, y=accuracy)) +
  geom_bar(stat="identity", position=position_dodge(), fill="#2a347a")+ 
  geom_text(aes(label=total), vjust=0.5, hjust=1.5, size=3.5, angle=90, color="white")+
  theme_classic()+
   theme(legend.position = "none", axis.text.x = element_text(size = 10))+
   scale_x_discrete(guide = guide_axis(angle = 90))+
  xlab("")+
  ggtitle('HGSVC HiFi')+
  ylab("Accuracy")

p3<- ggplot(data=df3, aes(x=gene, y=accuracy)) +
  geom_bar(stat="identity", position=position_dodge(), fill="#2a347a")+ 
  geom_text(aes(label=total), vjust=0.5, hjust=1.5, size=3.5, angle=90, color="white")+
  theme_classic()+
   theme(legend.position = "none", axis.text.x = element_text(size = 10))+
   scale_x_discrete(guide = guide_axis(angle = 90))+
  xlab("")+
  ggtitle('HPRC ONT')+
  ylab("Accuracy")

p4<- ggplot(data=df4, aes(x=gene, y=accuracy)) +
  geom_bar(stat="identity", position=position_dodge(), fill="#2a347a")+ 
  geom_text(aes(label=total), vjust=0.5, hjust=1.5, size=3.5, angle=90, color="white")+
  theme_classic()+
   theme(legend.position = "none", axis.text.x = element_text(size = 10))+
   scale_x_discrete(guide = guide_axis(angle = 90))+
  xlab("")+
  ggtitle('HGSVC CLR')+
  ylab("Accuracy")

# p3<- ggplot(data=df3, aes(x=gene, y=total)) +
#   geom_bar(stat="identity", position=position_dodge(), fill="#9fc3d5")+ #"#E69F00"
#   # theme_minimal()+
#   theme_classic()+
#    theme(legend.position = "none", axis.text.x = element_text(size = 10))+
#    scale_x_discrete(guide = guide_axis(angle = 90))+
#   xlab("")+
#    ggtitle('HPRC ONT')+
#   ylab("Count")+scale_fill_brewer(palette="Dark2")

# p4<- ggplot(data=df4, aes(x=gene, y=total)) +
#   geom_bar(stat="identity", position=position_dodge(), fill="#9fc3d5")+ #"#E69F00"
#   # theme_minimal()+
#   theme_classic()+
#    theme(legend.position = "none", axis.text.x = element_text(size = 10))+
#    scale_x_discrete(guide = guide_axis(angle = 90))+
#   xlab("")+
#  ggtitle('HGSVC CLR')+
# ylab("Count")+scale_fill_brewer(palette="Dark2")


prow <- plot_grid(
  p1 + theme(legend.position="none"),
  p2 + theme(legend.position="none"),
  p3 + theme(legend.position="none"),
  p4 + theme(legend.position="none"),
  align = 'vh',
  hjust = -1,
  nrow = 1
)
legend <- get_legend(
  p1 +
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom")
)
# legend = get_plot_component(p4, 'guide-box-top', return_all = TRUE)
plot_grid(prow, legend, ncol = 1, rel_heights = c(1, .1))
dev.off()