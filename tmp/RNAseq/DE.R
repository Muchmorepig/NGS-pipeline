library(ComplexHeatmap)
library(ggplot2)
library(DESeq2)
library(dplyr)

coldat <- read.csv("./dataSheet/RNAseq_Info.csv")
coldat <- coldat[!coldat$sample %in% c("F.S", "F.B", "M.S", "M.B"), ]

cts <- read.csv(
  "./upStream/timeSeries_RNA/featureCounts/tak12_v5_counts_noSB.csv",
  header = TRUE,
  row.names = 1
)
# 去掉低丰度基因
cts <- cts[rowSums(cts > 1) > 1, ]

# create DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = cts,
  colData = coldat,
  design = ~sample
)

# cor & PCA ----
corre <- cor(cts)

col_fun <- circlize::colorRamp2(
  c(min(corre), median(corre), max(corre)),
  c("#3bb0f0", "#ebebeb", "#f95426")
)

p1 <- Heatmap(
  corre,
  col = col_fun,
  name = "Correlation",
  show_column_names = FALSE,
  show_column_dend = FALSE,
  row_names_side = "left",
  heatmap_height = unit(4, "mm") * nrow(corre),
  heatmap_width = unit(4.4, "mm") * ncol(corre),
  rect_gp = gpar(col = "white")
)

source("./script/utils/calcHTsize_based_Complexheatmap.R")
(size <- Calc_HTSize(p1))

pdf("./downStrem/RNAres/plot/cor.pdf", width = size[1], height = size[2])
p1
dev.off()

vsd <- vst(dds)

df <- plotPCA(
  vsd,
  intgroup = c("sample"),
  returnData = TRUE,
  ntop = ifelse(0.1 * nrow(cts) > 3000, 0.1 * nrow(cts), 3000)
)
df$Gender <- stringr::str_split(df$group, "[.]", simplify = T)[, 1]
df$Stage <- stringr::str_split(df$group, "[.]", simplify = T)[, 2]

df2 <- df %>%
  group_by(group) %>%
  summarise(PC1 = mean(PC1), PC2 = mean(PC2))
df2$group <- stringr::str_split(df2$group, "[.]", simplify = T)[, 2]

p2 <- ggplot(df, aes(x = PC1, y = PC2)) +
  geom_point(aes(shape = Gender, color = Stage), stat = "unique", size = 4, alpha = 0.9) +
  ggrepel::geom_text_repel(
    data = df2,
    aes(label = group),
    size = 3.5,
    bg.r = 0.25,
    color = "#000103",
    bg.color = "#f5f0f0",
    seed = 1,
    box.padding = unit(0.5, "lines"),
    point.padding = unit(1, "lines")
  ) +
  labs(title = "PCA analysis", caption = "Map to V5") +
  theme_minimal() +
  coord_fixed()

ggsave(plot = p2, filename = "./downStrem/RNAres/plot/pca_noSB_notop2000.pdf", width = 8, height = 6)

# one-one res-----
source("./script/utils/makeCompareInfo.R")
comp <- CompInfo(colData = coldat, condition = "sample")
comp <- comp[c(6, 17, 27, 36, 44, 51), ]


dds <- DESeq(dds)

count_norm <- counts(dds, normalize = TRUE)
save(count_norm, coldat, file = "./downStrem/RNAres/dataSheet/count_norm_v5.rda")

source("./script/utils/one_one.R")

apply(comp, 1, one_one,
      Type = "sample",
      cutoff_padj = 0.05,
      od = "./downStrem/RNAres/diff",
      seq = "RNA"
)
