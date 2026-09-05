library(ggplot2)
library(dplyr)
library(patchwork)

g <- read.table("depth_G.txt", header = FALSE)
w <- read.table("depth_W.txt", header = FALSE)

colnames(g) <- c("chr", "pos", "depth")
colnames(w) <- c("chr", "pos", "depth")


g$sample <- "G"
w$sample <- "W"
df <- rbind(g, w)
df <- df %>% filter(pos %% 100 == 0)
df$zero_depth <- ifelse(df$depth == 0, TRUE, FALSE)


plot_depth <- function(data, sample_label, color_normal) {
  ggplot(data, aes(x = pos / 1e6, y = depth)) +
    geom_point(
      aes(color = zero_depth),
      size = 0.8,
      alpha = 0.8
    ) +
    scale_color_manual(
      values = c("FALSE" = color_normal, "TRUE" = "red"),
      labels = c("Depth > 0", "Depth = 0"),
      name = NULL
    ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "top",
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 11)
    ) +
    labs(
      x = "Chromosome D02 position (Mb)",
      y = "Depth",
      title = paste0(sample_label, "_chrD02")
    )
}


p1 <- plot_depth(df %>% filter(sample == "G"), "G", "#1F77B4")  # 蓝
p2 <- plot_depth(df %>% filter(sample == "W"), "W", "#2CA02C")  # 绿


p_final <- p1 / p2  # 上下排列
ggsave("Depth_D02_colored.png", p_final, width = 10, height = 6, dpi = 300)
