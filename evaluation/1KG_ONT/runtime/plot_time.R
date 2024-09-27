library(ggplot2)
library(grid)
library(cowplot)
library(dplyr)


pdf(file="time_dot.pdf", width=15, height=4, onefile=FALSE)


df<-read.table("kir_time.csv", sep=",", header=TRUE)
p1<-ggplot(df, aes(x=time, y=mem)) +
  geom_point()+
  xlab("CPU Time (h)") +
  ylab("Peak Memory (GB)")+
  ggtitle("KIR")
  # + 
  # geom_text(label=df$Allele,size=2)
## count mean value of time and mem, and wall_clock_time
print(paste("KIR mean time: ", mean(df$time), median(df$time)))
print(paste("mean mem: ", mean(df$mem), median(df$mem)))
print(paste("mean wall_clock_time: ", mean(df$wall_clock_time), median(df$wall_clock_time)))


df<-read.table("vdj_time.csv", sep=",", header=TRUE)
p2<-ggplot(df, aes(x=time, y=mem)) +
  geom_point()+
  xlab("CPU Time (h)") +
  ylab("Peak Memory (GB)")+
  ggtitle("IG+TCR")
  # + 
  # geom_text(label=df$Allele,size=2)
## count mean value of time and mem, and wall_clock_time
print(paste("IG+TCR mean time: ", mean(df$time), median(df$time)))
print(paste("mean mem: ", mean(df$mem), median(df$mem)))
print(paste("mean wall_clock_time: ", mean(df$wall_clock_time), median(df$wall_clock_time)))


df<-read.table("cyp_time.csv", sep=",", header=TRUE)
p3<-ggplot(df, aes(x=time, y=mem)) +
  geom_point()+
  xlab("CPU Time (h)") +
  ylab("Peak Memory (GB)")+
  ggtitle("CYP")
  # + 
  # geom_text(label=df$Allele,size=2)
## count mean value of time and mem, and wall_clock_time
print(paste("CYP mean time: ", mean(df$time), median(df$time)))
print(paste("mean mem: ", mean(df$mem), median(df$mem)))
print(paste("mean wall_clock_time: ", mean(df$wall_clock_time), median(df$wall_clock_time)))



prow <- plot_grid(
  p1 + theme(legend.position="none"),
  p2 + theme(legend.position="none"),
  p3 + theme(legend.position="none"),
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
