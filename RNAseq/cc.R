library(magrittr)
library(dplyr)
library(ComplexHeatmap)
source("~/script/tips_func/mf_Zscore.R")
source("~/script/tips_func/mf_mean.R")
mf_Minmax <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}



cell_cycle <- read.csv("~/reference/TAIR10/Arabidopsis_cell_cycle_genes_fx.csv")

cc_core <- read.table("/data5/wmc_data/Mpoly/cc_core_tair.txt", header = FALSE, sep = ",")
colnames(cc_core) <- c("sb1", "sb2", "TAIR")
anno <- readr::read_csv(
  "~/reference/Mpolymorpha/MpV5_tair_anno.csv",
  show_col_types = FALSE
)


mp_cc <- dplyr::select(anno, geneId, TAIR) %>% inner_join(cell_cycle)
mp_cc <- dplyr::select(anno, geneId, TAIR, symbol) %>% inner_join(cc_core)
mp_cc <- unique(mp_cc)

counts <- read.csv("./result/count_norm.csv", row.names = 1)
head(counts)

(colData <- data.frame(
  name = colnames(counts),
  group = as.factor(rep(c("C1d", "C3d", "W1d", "W3d"), c(2,2,3,3)))
))

counts <- counts[mp_cc$geneId, ]

dat <- mf_mean1(counts, colData, "group") %>% mf_zscore()

col_fun <- circlize::colorRamp2(c(-2, 0, 2), c("#2fa1dd", "white", "#f87669"))

tmp <- as.data.frame(mp_cc)[mp_cc$geneId %in% rownames(dat), ]

p <- Heatmap(dat,
  name = "Z-Score\nAve Expr",
  col = col_fun,
  # row_split = Group,
  row_title = NULL,
  row_gap = unit(3, "mm"),
  rect_gp = gpar(col = "gray", lwd = 1.5),
  row_names_side = "left",
  column_names_rot = 0,
  column_names_centered = TRUE,
  top_annotation = HeatmapAnnotation(Sum_MinMax = anno_barplot(mf_Minmax(rowSums(t(dat))))),
  right_annotation = rowAnnotation(
    foo = anno_text(tmp$symbol)
  )
)



pdf("./result/plot/cc.pdf", width = 5,height = 11)
p
dev.off()

# counts <- read.csv("./result/count_norm.csv", row.names = 1)
# head(counts)
# 
# dat2 <- edgeR::cpm(counts)
# dat2 <- log10(dat2 + 1)
# dat2 <- dat2[unique(mp_cc$geneId), ]
# dat2 <- mf_mean1(dat2, colData, "group")
# dat2 <- dat2[rowSums(dat2) > 0, ]
# 
# 
# tmp2 <- as.data.frame(mp_cc)[unique(mp_cc$geneId) %in% rownames(dat2), ]
# 
# col_fun2 <- circlize::colorRamp2(c(0, 2), c("white", "#f87669"))
# 
# Heatmap(dat2,
#         name = "Z-Score\nAve Expr",
#         col = col_fun2,
#         # row_split = Group,
#         row_title = NULL,
#         row_gap = unit(3, "mm"),
#         rect_gp = gpar(col = "gray", lwd = 1.5),
#         row_names_side = "left",
#         column_names_rot = 0,
#         column_names_centered = TRUE,
#         top_annotation = HeatmapAnnotation(Sum_MinMax = anno_barplot(mf_Minmax(rowSums(t(dat2))))),
#         right_annotation = rowAnnotation(
#           foo = anno_text(tmp2$symbol)
#         )
# )
