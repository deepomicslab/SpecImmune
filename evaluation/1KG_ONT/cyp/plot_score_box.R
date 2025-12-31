## load the csv file as a matrix
library(ggplot2)
library(cowplot)
library(dplyr)
library(ggpubr)


df<-read.table("cyp_activity.csv", sep=",", header=TRUE)

# Add sample size to x-axis labels
df <- df %>%
  group_by(Super_Pop) %>%
  mutate(Super_Pop_label = paste0(Super_Pop, " (n=", n(), ")")) %>%
  ungroup()

## EAS vs non-EAS comparison
df1<-filter(df, Super_Pop == "EAS")
df2<-filter(df, Super_Pop != "EAS")
t_result <- t.test(df1$Activity_score, df2$Activity_score)
print(t_result)

# Format p-value for display with detailed precision
p_val <- t_result$p.value
p_label <- sprintf("t-test, p = %.2e", p_val)

# Create comparison dataframe
df_comparison <- df %>%
  mutate(Group = ifelse(Super_Pop == "EAS", "EAS", "Non-EAS"))

# Figure 1: All super populations
pdf(file="cyp_activity_box.pdf", width=4, height=4, onefile=FALSE)

p1 <- ggplot(df, aes(x=Super_Pop_label, y = Activity_score, fill=Super_Pop)) +
  geom_boxplot(outlier.shape = NA, alpha=0.7)  +
  geom_jitter(width=0.2, size=0.8, alpha=0.4) +  # Show individual data points
  stat_summary(fun=mean, geom="point", shape=20, size=3, color="skyblue", fill="skyblue") +
  theme_classic()+
  xlab('')+
  ylab('CYP2D6 Activity Score')+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 9)) +
  scale_fill_manual(values = c("#d9e6eb", "#9fc3d5", "#8f96bd", "#2a347a", "#d6d69b"))

print(p1)

dev.off()

# Figure 2: EAS vs Non-EAS with p-value
pdf(file="cyp_activity_box_eas_vs_noneas.pdf", width=2, height=2.5, onefile=FALSE)

p2 <- ggplot(df_comparison, aes(x=Group, y = Activity_score, fill=Group)) +
  geom_boxplot(outlier.shape = NA, alpha=0.7)  +
  geom_jitter(width=0.2, size=0.8, alpha=0.4) +
  stat_summary(fun=mean, geom="point", shape=20, size=3, color="skyblue", fill="skyblue") +
  theme_classic()+
  xlab('')+
  ylab('CYP2D6 Activity Score')+
  theme(legend.position = "none") +
  scale_fill_manual(values = c("#8f96bd", "#d6d69b")) +
  stat_compare_means(method = "t.test", label.y = 3.5, size = 3)

print(p2)

dev.off()