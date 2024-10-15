library(ggplot2)

df <- data.frame(method=c("pangu", "SpecLong"),
                 Accuracy=c(0, 0.95))

pdf("figures/cyp_amplicon.pdf", width=2.2, height=4)
# df$Group <- factor(df$Group, levels=c("Thomas-Abun", "Hybrid"))


 p <- ggplot(data=df, aes(x=method, y=Accuracy, fill=method))+
   geom_bar(stat="identity", color="black", position=position_dodge(), width=0.6)+
   theme_classic()+
   scale_x_discrete(guide = guide_axis(angle = 30))+
   coord_cartesian(ylim=c(0,1))+
   xlab("")+
   theme(legend.position = "none")+
  #  theme(legend.position = "none", axis.text.x = element_text(size = 13))+
ggtitle("PacBio Amplicon")
#  p + scale_fill_manual(values=c("slateblue1" , "tomato"))
#  p + scale_fill_manual(values=c("#d9e6eb" , "#9fc3d5"))
  p + scale_fill_manual(values=c("#9fc3d5","#2a347a"))
dev.off()