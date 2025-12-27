library(ggplot2)
library(grid)
library(cowplot)
library(dplyr)

pdf(file="figures/cyp_depth.pdf", width=10, height=3, onefile=FALSE)

### -------- 1kGP --------
df <- read.table("cyp_depth_cutoff.csv", sep=",", header=TRUE)

# Build axis labels (one N per cutoff)
axis_df <- df %>%
  group_by(cutoff) %>%
  summarize(N = max(total), .groups="drop")

p1 <- ggplot(data=df, aes(x=cutoff, y=accuracy, color=method)) +
  geom_line() +
  geom_point() +
  xlab("Minimum reads per sample (filter cutoff)") +
  ylab("Accuracy") +
  ggtitle("1kGP") +
  scale_x_continuous(
    breaks = axis_df$cutoff,
    labels = axis_df$cutoff,
    sec.axis = sec_axis(
      trans = ~ .,
      breaks = axis_df$cutoff,
      labels = axis_df$N,
      name = "Alleles retained (N)"
    )
  ) +
  scale_color_manual(values = c("#9fc3d5","#2a347a")) +
  theme_classic()

### -------- HPRC HiFi --------
df <- read.table("hprc_hifi_cyp_depth_cutoff.csv", sep=",", header=TRUE)

axis_df <- df %>%
  group_by(cutoff) %>%
  summarize(N = max(total), .groups="drop")

p2 <- ggplot(data=df, aes(x=cutoff, y=accuracy, color=method)) +
  geom_line() +
  geom_point() +
  xlab("Minimum reads per sample (filter cutoff)") +
  ylab("Accuracy") +
  ggtitle("HPRC HiFi") +
  scale_x_continuous(
    breaks = axis_df$cutoff,
    labels = axis_df$cutoff,
    sec.axis = sec_axis(
      trans = ~ .,
      breaks = axis_df$cutoff,
      labels = axis_df$N,
      name = "Alleles retained (N)"
    )
  ) +
  scale_color_manual(values = c("#9fc3d5","#2a347a")) +
  theme_classic()

### -------- HPRC ONT --------
df <- read.table("hprc_ont_cyp_depth_cutoff.csv", sep=",", header=TRUE)

axis_df <- df %>%
  group_by(cutoff) %>%
  summarize(N = max(total), .groups="drop")

p3 <- ggplot(data=df, aes(x=cutoff, y=accuracy, color=method)) +
  geom_line() +
  geom_point() +
  xlab("Minimum reads per sample (filter cutoff)") +
  ylab("Accuracy") +
  ggtitle("HPRC ONT") +
  scale_x_continuous(
    breaks = axis_df$cutoff,
    labels = axis_df$cutoff,
    sec.axis = sec_axis(
      trans = ~ .,
      breaks = axis_df$cutoff,
      labels = axis_df$N,
      name = "Alleles retained (N)"
    )
  ) +
  scale_color_manual(values = c("#9fc3d5","#2a347a")) +
  theme_classic()

### -------- Combine panels --------
prow <- plot_grid(
  p2 + theme(legend.position="none"),
  p3 + theme(legend.position="none"),
  p1 + theme(legend.position="none"),
  align = 'vh',
  hjust = -1,
  nrow = 1
)

legend <- get_legend(
  p1 +
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom")
)

plot_grid(prow, legend, ncol = 1, rel_heights = c(1, .1))

dev.off()
