library(ggplot2)
library(ggpubr)
library(grid)
library(cowplot)
library(dplyr)

df <- read.table("vdj_snp_count.csv", sep=",", header=TRUE)
pdf(file="vdj_snp.pdf", width=13, height=7, onefile=FALSE)


df1<-filter(df,Group=="TRA")
p1 <- ggplot(df1, aes(x=new_pop,y=Hete_Variant_Num,fill=new_pop)) + 
        geom_boxplot() + 
        # geom_bar(stat='summary') + 
    stat_summary(fun=mean, geom="point", shape=20, size=5, color="yellow", fill="yellow") +
        theme_bw()+
        theme(axis.ticks=element_blank(), axis.title=element_blank(),axis.text.x = element_blank(), panel.grid  = element_blank())+
         scale_fill_brewer(palette="Dark2")+ stat_compare_means(method = "t.test")+
         ggtitle("TRA")

df1<-filter(df,Group=="TRD")
p2 <- ggplot(df1, aes(x=new_pop,y=Hete_Variant_Num,fill=new_pop)) + 
        geom_boxplot() + 
        theme_bw()+
        theme(axis.ticks=element_blank(), axis.title=element_blank(),axis.text.x = element_blank(), panel.grid  = element_blank())+
         scale_fill_brewer(palette="Dark2")+ stat_compare_means(method = "t.test")+
         ggtitle("TRD")

df1<-filter(df,Group=="TRB")
p3 <- ggplot(df1, aes(x=new_pop,y=Hete_Variant_Num,fill=new_pop)) + 
        geom_boxplot() + 
        theme_bw()+
        theme(axis.ticks=element_blank(), axis.title=element_blank(),axis.text.x = element_blank(), panel.grid  = element_blank())+
         scale_fill_brewer(palette="Dark2")+ stat_compare_means(method = "t.test")+
         ggtitle("TRB")

df1<-filter(df,Group=="TRG")
p4 <- ggplot(df1, aes(x=new_pop,y=Hete_Variant_Num,fill=new_pop)) + 
        geom_boxplot() + 
        theme_bw()+
        theme(axis.ticks=element_blank(), axis.title=element_blank(),axis.text.x = element_blank(), panel.grid  = element_blank())+
         scale_fill_brewer(palette="Dark2")+ stat_compare_means(method = "t.test")+
         ggtitle("TRG")


df1<-filter(df,Group=="IGH")
p5 <- ggplot(df1, aes(x=new_pop,y=Hete_Variant_Num,fill=new_pop)) + 
        geom_boxplot() + 
        theme_bw()+
        theme(axis.ticks=element_blank(), axis.title=element_blank(),axis.text.x = element_blank(), panel.grid  = element_blank())+
         scale_fill_brewer(palette="Dark2")+ stat_compare_means(method = "t.test")+
         ggtitle("IGH")


df1<-filter(df,Group=="IGK")
p6 <- ggplot(df1, aes(x=new_pop,y=Hete_Variant_Num,fill=new_pop)) + 
        geom_boxplot() + 
        theme_bw()+
        theme(axis.ticks=element_blank(), axis.title=element_blank(),axis.text.x = element_blank(), panel.grid  = element_blank())+
         scale_fill_brewer(palette="Dark2")+ stat_compare_means(method = "t.test")+
         ggtitle("IGK")


df1<-filter(df,Group=="IGL")
p7 <- ggplot(df1, aes(x=new_pop,y=Hete_Variant_Num,fill=new_pop)) + 
        geom_boxplot() + 
        theme_bw()+
        theme(axis.ticks=element_blank(), axis.title=element_blank(),axis.text.x = element_blank(), panel.grid  = element_blank())+
         scale_fill_brewer(palette="Dark2")+ stat_compare_means(method = "t.test")+
         ggtitle("IGL")


  # p1
prow <- plot_grid(
  p1 + theme(legend.position="none"),
  p2 + theme(legend.position="none"),
  p3 + theme(legend.position="none"),
  p4 + theme(legend.position="none"),
  p5 + theme(legend.position="none"),
  p6 + theme(legend.position="none"),
  p7 + theme(legend.position="none"),
  align = 'vh',
  hjust = -1,
  nrow = 2,
  ncol = 4
)

legend <- get_legend(
  p1 +
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom")
)

plot_grid(prow, legend, ncol = 1, rel_heights = c(1, .1))



dev.off()