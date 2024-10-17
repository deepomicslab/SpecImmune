library(ggplot2)
library(ggpubr)
library(ggplot2)
library(grid)
library(cowplot)

pdf(file="figures/hla_class_kir_vdj.pdf", width=8.2, height=6, onefile=FALSE)
df<-read.table("tree/vdj_cyp_hla_kir.csv", sep=",", header=TRUE)
## sort the df by ALD in reverse order
# df <- df[order(df$ALD, decreasing = TRUE),]
## maintain the order of the type
## output all types
unique(df$type)

df$type <- factor(df$type, levels = c("HLA", "KIR", "IG", "TCR", "(HLA,KIR)", "(CYP,HLA)", "(CYP,KIR)", "(IG,TCR)", "(IG,KIR)", "(HLA,IG)", "(HLA,TCR)", "(CYP,IG)", "(CYP,TCR)","(KIR,TCR)"))
# "TCR"       "IG"        "HLA"       "KIR"       "(IG,TCR)"  "(HLA,KIR)"
# "(CYP,KIR)" "(CYP,HLA)" "(HLA,TCR)" "(HLA,IG)"  "(IG,KIR)"  "(CYP,IG)" 
# "(KIR,TCR)" "(CYP,TCR)"

  vysg<-ggplot(df, aes(x=ALD, y = after_stat(density))) +
    # geom_histogram(position="identity", alpha=0.5,bins = 50, fill="skyblue", color="skyblue")+
    geom_histogram(aes(y=after_stat(density)), 
                   position="identity", color="#9fc3d5", fill="#9fc3d5", bins=100)+
    xlab("Wmin")+
    ylab("Count")+
    # geom_density(lwd = 0.3, colour = "#F8766D",
    #              fill = 4, alpha = 0.25)+
    geom_density(alpha=.2, color="#2a347a") +
    # xlim(0, 2000)+
    theme_bw()

vysg<-vysg+facet_wrap(~ type , ncol = 4, scales = "free")+
xlab('minALD')+
  theme(legend.position = "none")

## select the df with type=="HLA"
df2 <- df[df$type == "(HLA,IG)",]
## count how many loci in df2
nrow(df2)
## count how many loci with ALD > 0.2 in df2
nrow(df2[df2$ALD < 0.2,])
### count proportion of loci with ALD > 0.2 in df2
nrow(df2[df2$ALD < 0.2,])/nrow(df2)
## sort df2 by ALD in reverse order, and head 10
df2 <- df2[order(df2$ALD, decreasing = TRUE),]
head(df2, 10)





vysg
