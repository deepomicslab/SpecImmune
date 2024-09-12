library(ggplot2)
library(ggpubr)
library(grid)
library(cowplot)
library(dplyr)

df <- read.table("vdj_snp_count.csv", sep=",", header=TRUE)
pdf(file="vdj_snp.pdf", width=20, height=120, onefile=FALSE)






vysg <- ggplot(df, aes(x=new_pop,y=Hete_Variant_Num,fill=new_pop)) + 
          geom_boxplot() + 
        theme_bw()+
        theme(axis.ticks=element_blank(), axis.title=element_blank(),axis.text.x = element_blank(), panel.grid  = element_blank())+
         scale_fill_brewer(palette="Dark2")+ stat_compare_means(method = "t.test")


# vysg <- ggplot(df, aes(x=new_pop,y=Hete_Variant_Num,fill=new_pop)) + 
#           # geom_boxplot() + 
#           geom_bar(stat='summary') + 
#         theme_bw()+
#         theme(axis.ticks=element_blank(), axis.title=element_blank(),axis.text.x = element_blank(), panel.grid  = element_blank())+
#          scale_fill_brewer(palette="Dark2")+ stat_compare_means(method = "t.test")+
#          ylim(c(0,30))

vysg<-vysg+facet_wrap(~ Gene, ncol=10, scales = "free")
vysg

dev.off()