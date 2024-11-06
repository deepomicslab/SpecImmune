library(ggplot2)
# library(grid)
# library(cowplot)
# library(dplyr)


pdf(file="figures/cyp_sim.pdf", width=4, height=3, onefile=FALSE)


df<-read.table("sim_all_gene.csv", sep=",", header=TRUE)
df
## in relace the specLong to spec in the df
# df$method <- gsub("SpecLong", "SpecImmune", df$method)
p1<-ggplot(data=df, aes(x=Gene, y=Accuracy)) +
  geom_bar(stat="identity", color="#9fc3d5", fill="#9fc3d5")+
  xlab("")+
  ylab("Accuracy")+
  scale_color_manual(values = c("#9fc3d5","#2a347a"))+theme_classic()+
   scale_x_discrete(guide = guide_axis(angle = 90))
p1
dev.off()