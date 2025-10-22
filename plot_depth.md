# 1 对 D02 染色体计算深度
```bash
samtools depth -r D02 1077bulk-G.final.bam > depth_G.txt
samtools depth -r D02 1077bulk-W.final.bam > depth_W.txt
```

# 2 R可视化
## 2.1 绘制整条染色体
```R
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

```

## 2.2 绘制特定染色体
```R
library(ggplot2)
library(dplyr)
library(patchwork)


region_start <- 62965000
region_end <- 69769000


g <- read.table("depth_G.txt", header = FALSE)
w <- read.table("depth_W.txt", header = FALSE)
colnames(g) <- c("chr", "pos", "depth")
colnames(w) <- c("chr", "pos", "depth")
g$sample <- "G"
w$sample <- "W"

df <- rbind(g, w)


df <- df %>%
  filter(pos >= region_start & pos <= region_end)


df <- df %>% filter(pos %% 100 == 0)


plot_depth <- function(data, sample_label, color_normal) {
  ggplot(data, aes(x = pos / 1e3, y = depth)) +   # 改为kb
    geom_rect(
      aes(xmin = region_start / 1e3, xmax = region_end / 1e3,
          ymin = -Inf, ymax = Inf),
      fill = "grey85", alpha = 0.5, inherit.aes = FALSE
    ) +
    geom_point(
      data = subset(data, depth > 0),
      color = color_normal,
      size = 0.7,
      alpha = 0.6
    ) +
    geom_point(
      data = subset(data, depth == 0),
      color = "red",
      size = 1.2,
      alpha = 0.9
    ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 11)
    ) +
    labs(
      x = "Chromosome D02 position (Kb)",      # 修改坐标标签
      y = "Depth",
      title = paste0(sample_label, "_chrD02 (", region_start, "-", region_end, ")")
    )
}


p1 <- plot_depth(df %>% filter(sample == "G"), "G", "#1F77B4")
p2 <- plot_depth(df %>% filter(sample == "W"), "W", "#2CA02C")

p_final <- p1 / p2

ggsave("Depth_D02_region_highlight.png", p_final, width = 10, height = 6, dpi = 300)
```
