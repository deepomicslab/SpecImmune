library(ggplot2)
library(grid)
library(cowplot)
library(dplyr)

# df<-read.table("HPRC_ont_read_num.csv", sep=",", header=TRUE)   #12 / 37
# df<-read.table("HGSCV2_hifi_read_num.csv", sep=",", header=TRUE)  #9/9
df<-read.table("hgscv2_clr_read_num.csv", sep=",", header=TRUE)  # 16/19
# df<-read.table("HPRC_hifi_read_num.csv", sep=",", header=TRUE)  #40/45
df <- dplyr::filter(df, gene == "KIR3DP1")
## count the number of elements in df
nrow(df)
df <- dplyr::filter(df, read_num >=5)
nrow(df)