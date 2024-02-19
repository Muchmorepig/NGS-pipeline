library(dplyr)
library(tidyr)
library(cluster)
library(ggplot2)
library(ComplexHeatmap)

# selece Peaks for k-mean-------------
(diff_files <- list.files("./downStrem/RNAres/diff", full.names = TRUE))

diff_list <- lapply(diff_files, function(x) {
  read.csv(x) %>%
    subset(padj < 0.05 & abs(log2FoldChange) > 1) %>%
    mutate(change = ifelse(log2FoldChange > 1, "up", "down"))
})

tmp <- do.call(rbind, diff_list) %>% select(geneId, VS, change)

tmp2 <- tmp %>%
  group_by(VS) %>%
  nest() %>%
  rowwise() %>%
  mutate(tb = list(table(data$change)))

names(tmp2$tb) <- tmp2$VS

pp <- do.call(cbind, tmp2$tb) %>%
  as_tibble(rownames = "Change") %>%
  pivot_longer(., 2:7, names_to = "VS") %>%
  mutate(
    x = factor(VS,
      levels = c(
        "F.0d_VS_M.0d", "F.1d_VS_M.1d", "F.3d_VS_M.3d",
        "F.6d_VS_M.6d", "F.11d_VS_M.11d", "F.20d_VS_M.20d"
      )
    ),
    y = ifelse(Change == "down", -value, value)
  ) %>%
  ggplot(aes(x = x, y = y, fill = Change)) +
  geom_bar(position = "stack", stat = "identity") +
  theme_light() +
  theme(
    axis.title = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.title = element_blank()
  ) +
  scale_fill_manual(
    breaks = c("down", "up"),
    values = c("#2fa1dd", "#f87669")
  ) +
  scale_y_continuous(
    breaks = c(600, 300, 0, -300, -600),
    labels = c("600", "300", "0", "300", "600")
  ) +
  scale_x_discrete(labels = paste0(c(0, 1, 3, 6, 11, 20), "d")) +
  labs(caption = "Female vs Male, LFC's absolute value over 1")

ggsave(pp, filename = "./downStrem/RNAres/plot/diff_stat.pdf", width = 8, height = 6)

sel_gene <- unique(tmp$geneId)


load("./downStrem/RNAres/dataSheet/count_norm_v5.rda")

# peak_norm <- log1p(peak_norm)

count_sel <- count_norm[sel_gene, ]

source("~/script/tips_func/mf_mean.R")
source("~/script/tips_func/mf_Zscore.R")

count_sel_mean <- mf_mean2(count_sel, group = gsub("\\.[1-3]$", "", colnames(count_sel)))
count_sel_mean_zscore <- mf_zscore(count_sel_mean)

# kmean ---------
set.seed(1314)

km <- kmeans(count_sel_mean_zscore, centers = 12, iter.max = 50)
km_res <- cbind(cluster = km$cluster, count_sel_mean_zscore)
pd <- km_res

save(km, file = "./downStrem/RNAres/dataSheet/sel_km.rda")
write.csv(pd, "./downStrem/RNAres/dataSheet/sel_km_zscore.csv")

col_fun <- circlize::colorRamp2(
  c(min(pd[, -1]), 0, max(pd[, -1])),
  c("#0080ff", "#e8e8e8", "#ee0038")
)

Group <- factor(rep(c("Tak1", "Tak2"), each = 6), levels = c("Tak1", "Tak2"))

top_annotation <- HeatmapAnnotation(
  cluster = anno_block(
    gp = gpar(fill = c("#2fa1dd", "#f87669")),
    labels = c("Tak1", "Tak2"),
    labels_gp = gpar(col = "white", fontsize = 12)
  )
)

h1 <- Heatmap(
  pd[, -1],
  col = col_fun,
  top_annotation = top_annotation,
  column_split = Group,
  column_title = NULL,
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  show_row_names = FALSE,
  row_split = factor(pd[, "cluster"], levels = sort(unique(pd[, "cluster"]))),
  use_raster = TRUE,
  heatmap_legend_param = list(
    at = c(round(min(pd[, -1])), 0, round(max(pd[, -1]))),
    direction = "horizontal",
    title_position = "leftcenter",
    legend_width = unit(3.8, "cm"),
    title = "Z-Score"
  )
)

pdf("./downStrem/RNAres/plot/km_hp.pdf", height = 12, width = 4)
ComplexHeatmap::draw(h1, heatmap_legend_side = "top")
dev.off()

pd2 <- pd %>%
  as_tibble(rownames = "geneId") %>%
  pivot_longer(3:14, names_to = "stage", values_to = "zscore") %>%
  mutate(gender = gsub("\\..+d", "", .$stage), day = gsub("[FM]\\.", "", .$stage))

pd2$cluster <- factor(pd2$cluster)
pd2$x <- rep(c(1:6), times = length(sel_gene) * 2)


.mytheme <- theme_minimal() + theme(
  # axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  # axis.ticks = element_blank(),
  axis.title = element_blank(),
  # strip.text = element_blank(),
  # panel.spacing.x = unit(-1, "lines"),
  panel.spacing.y = unit(0, "lines"),
  # panel.border = element_blank(),
  # axis.ticks = element_blank(),
  panel.grid.minor = element_blank(),
  panel.grid.major.x = element_blank(),
  legend.title = element_blank(),
  # strip.background = element_rect(color = "#f3eaea"),
  strip.text = element_text(face = "bold")
  # legend.position = "top",
  # legend.background = element_rect(
  #   color = "#f3eaea"
  # )
)

p1 <- ggplot(data = pd2, aes(x = x, y = zscore, color = gender)) +
  # geom_point() +
  stat_smooth(method = loess, linewidth = 1.2, span = .6, aes(color = gender)) +
  facet_grid(cluster ~ .) +
  scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 6), labels = c("0d", "1d", "3d", "6d", "11d", "20d")) +
  scale_color_manual(values = c("M" = "#2fa1dd", "F" = "#f87669")) +
  .mytheme +
  guides(color = "none")

pdf("./downStrem/RNAres/plot/km_line.pdf", height = 11.5, width = 1.8)
p1
dev.off()

alias <- read.csv("/data5/wmc_data/reference/Mpolymorpha/MpolyV5_Tair_anno.csv")

pd %>%
  as_tibble(rownames = "geneId") %>%
  select(geneId, cluster) %>%
  left_join(alias, multiple = "all") %>%
  write.csv("./downStrem/RNAres/dataSheet/sel_km_anno.csv")