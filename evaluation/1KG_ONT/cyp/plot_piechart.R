library(ggplot2)

df <- read.table("cyp_pheno.csv", sep=",", header=TRUE)
pdf(file="cyp_pheno.pdf", width=8, height=4, onefile=FALSE)
vysg <- ggplot(df, aes(x=1,y=Frequency,fill=Phenotype)) + 
          geom_bar(stat="identity",width=2) + 
          coord_polar(theta='y')+
        theme_bw()+
        theme(axis.ticks=element_blank(), axis.title=element_blank(),axis.text.y = element_blank(),axis.text.x = element_blank(), panel.grid  = element_blank())+
         scale_fill_brewer(palette="Dark2")
vysg<-vysg+facet_wrap(~ Super_Pop)
vysg

dev.off()