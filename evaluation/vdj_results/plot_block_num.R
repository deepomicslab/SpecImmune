  library(ggplot2)
#   library(dplyr)
#   require(gridExtra)
#   library(grid)
  library(cowplot)
  
  df<-read.table("block_loci_num.csv", sep=",", header=TRUE)
  pdf(file="figures/block_num.pdf", width=4, height=3, onefile=FALSE)
  
  p1<-ggplot(df, aes(x=Num, y = after_stat(density))) +
    # geom_histogram(position="identity", alpha=0.5,bins = 50, fill="skyblue", color="skyblue")+
    geom_histogram(aes(y=after_stat(density)), 
                   position="identity", color="blue", fill="lightblue")+
    xlab("No. of gene loci")+
    ylab("Density")+
    # geom_density(lwd = 0.3, colour = "#F8766D",
    #              fill = 4, alpha = 0.25)+
    geom_density(alpha=.2) +
    # xlim(0, 2000)+
    theme_classic()
#     +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.line = element_line(colour = "black"),legend.position = "top")+ geom_vline(aes(xintercept=mean(Bkp_count)),
                                                                                        # color="black", linetype="dashed", linewidth=1)
  
  p1
  
  dev.off()
  
  
