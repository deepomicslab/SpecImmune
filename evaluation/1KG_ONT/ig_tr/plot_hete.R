library(ggplot2)
library(ggpubr)
library(grid)
library(cowplot)
library(dplyr)

df <- read.table("vdj_hete_freq.csv", sep=",", header=TRUE)
pdf(file="vdj_hete.pdf", width=14, height=6, onefile=FALSE)


df1<-filter(df, Gene_Group2 == "IG")
p1 <- ggplot(df1, aes(x=Gene,y=Hete_freq,fill=Gene_Group)) + 
          geom_bar(stat="identity", width=0.5) +
  theme_classic()+
   theme(axis.text.x = element_text(size = 4))+
   scale_x_discrete(guide = guide_axis(angle = 90))+
  ylab("Heterozygous Frequency")+scale_fill_brewer(palette="Dark2")+xlab('')


df1<-filter(df, Gene_Group2 == "TR")
p2 <- ggplot(df1, aes(x=Gene,y=Hete_freq,fill=Gene_Group)) + 
          geom_bar(stat="identity", width=0.5) +
  theme_classic()+
   theme(axis.text.x = element_text(size = 4))+
   scale_x_discrete(guide = guide_axis(angle = 90))+
  ylab("Heterozygous Frequency")+scale_fill_brewer(palette="Dark2")+xlab('')

prow <- plot_grid(
  p1 + theme(legend.position="top"),
  p2 + theme(legend.position="top"),
  align = 'vh',
  hjust = -1,
  nrow = 2
)

plot_grid(prow, ncol = 1, rel_heights = c(1, .1))

# vysg <- ggplot(df, aes(x=Gene,y=Hete_freq)) + 
#           geom_bar(stat="identity") + 
  # vysg<-ggplot(df, aes(x=Hete_freq, y = after_stat(density))) +
  #   geom_histogram(aes(y=after_stat(density)), 
  #                  position="identity", color="white", fill="skyblue", alpha=0.5,bins = 20)+
  #       theme_bw()+
  #       theme(axis.ticks=element_blank(), axis.title=element_blank(), panel.grid  = element_blank())+
  #        scale_fill_brewer(palette="Dark2")



# vysg <- ggplot(df, aes(x=Gene,y=Hete_freq)) + 
#           geom_bar(stat="identity") + 
#         theme_bw()+
#         theme(axis.ticks=element_blank(), axis.title=element_blank(), panel.grid  = element_blank())+
#          scale_fill_brewer(palette="Dark2")


# vysg<-vysg+facet_wrap(~ Gene_Group, nrow = 2, scales = "free")
# vysg

dev.off()