library(ggplot2)
library(grid)
library(cowplot)
library(ggrepel)

pdf(file="figures/all_depth.pdf", width=5.5, height=4, onefile=FALSE)
df<-read.table("all_loci_depth.csv", sep=",", header=TRUE)

# Set factor levels for dataset to control legend order
df$dataset <- factor(df$dataset, levels = c("HPRC HiFi", "HPRC ONT", "HGSVC HiFi", "HGSVC CLR"))

# Add label combining accuracy and total alleles
df$label <- paste0(sprintf("%.1f%%", df$accuracy * 100), "\n(n=", df$total, ")")

ggplot(data=df, aes(x=depth, y=accuracy, color=dataset)) +
  geom_line(size=1)+
  geom_point()+
  geom_text_repel(aes(label=total), size=2.5, show.legend=FALSE, 
                  max.overlaps = 50, 
                  box.padding = 1, 
                  point.padding = 0.5,
                  min.segment.length = 0,
                  force = 5,
                  force_pull = 0.5,
                  seed = 42,
                  direction = "both",
                  segment.size = 0.3,
                  segment.alpha = 0.6) +  
  xlab("Minimum depth (filter cutoff)")+
  ylim(c(0.95,1.005))+  # More space for labels
  scale_x_continuous(breaks = c(0, 5, 10, 15, 20)) +
  geom_line(aes(linetype=dataset))+
  ylab("Accuracy")+theme_classic()+ scale_color_manual(values = c("#d9e6eb", "#9fc3d5", "#8f96bd", "#2a347a", "#d6d69b"))+
  theme(legend.position="top",
        plot.margin = margin(10, 10, 10, 10))
dev.off()
