## load the csv file as a matrix
library(ggplot2)
library(aplot)


pdf(file="hap_freq_heatmap.pdf", width=12, height=5, onefile=FALSE)


df<-read.table("cyp_hap_freq.csv", sep=",", header=TRUE)

## assign x order
df$Pop <- factor(df$Pop, levels = c('CDX', 'CHB', 'CHS', 'JPT', 'KHV', 'BEB', 'GIH', 'ITU', 'PJL', 'STU', 'ACB', 'ASW', 'ESN', 'GWD', 'LWK', 'MSL', 'YRI', 'CEU', 'FIN', 'GBR', 'IBS', 'TSI', 'CLM', 'MXL', 'PEL', 'PUR'))
# df$Allele<- factor(df$Allele, levels =c('*1', '*2', '*10', '*4', '*41', '*17', '*36+*10', '*35', '*141', '*29', '*171', '*5', '*68+*4', '*3', '*45', '*9', '*2x2', '*39', '*21', '*6', '*108', '*4x2', '*36', '*56', '*27', '*43', '*121', '*101', '*1x2', '*106', '*32', '*150', '*31', '*146', '*46', '*136', '*14', '*117', '*71', '*127', '*68+*2', '*82', '*129', '*7', '*28', '*68', '*4.013+*4', '*36x2+*10', '*112', '*120', '*38', '*111', '*59', '*132', '*119', '*83', '*10+*10', '*84', '*58', '*15', '*86', '*92', '*170', '*99', '*115', '*49', '*110', '*122', '*11', '*43x3', '*151', '*2x3', '*22', '*171x2', '*19', '*159', '*18', '*144', '*64', '*165', '*90', '*33', '*171x3', '*89', '*88', '*169', '*60', '*75', '*113', '*100', '*166', '*74', '*145', '*4.013', '*17x2', '*125'))



p1<-ggplot(data = df, aes(x=Allele, y=Pop, fill=Frequency)) + 
  geom_tile()+
   theme(axis.text.x = element_text(size = 8))+
   scale_x_discrete(guide = guide_axis(angle = 90))+
 geom_tile(color = "white")+
 scale_fill_gradient2(low = "#8fb4be", high = "#d93f49", mid = "white",midpoint = 0.1,
    space = "Lab",   #midpoint = 0.3, limit = c(0,0.5),
    name="Freq") +
#   theme_classic()

    theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank())

df2<-read.csv("../hla/hla_color.csv",header=T, sep=",")
df2$y<-factor(df2$y,levels = rev(df2$y))
p2<-ggplot(df2,aes(x=x,y=y))+
  geom_tile(aes(fill=group))+
  scale_x_continuous(expand = c(0,0))+
  theme(panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = "left",
        legend.title = element_blank())+
  scale_fill_manual(values = c("#d9e6eb", "#9fc3d5", "#8f96bd", "#2a347a", "#d6d69b"))
#         +
#   scale_fill_manual(values = c("green","blue","red"))
df3<-read.csv("cyp_hap_count.csv",header=T, sep=",")
df3$y<-factor(df3$Hap,levels = rev(df3$Hap))
p3<-ggplot(df3,aes(x=Hap,y=Count))+
  geom_bar(stat="identity")+
  theme(panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "left",
        legend.title = element_blank())

p <- p1%>%
  insert_left(p2,width = 0.02) %>%
  insert_top(p3, height = 0.2)
p

dev.off()