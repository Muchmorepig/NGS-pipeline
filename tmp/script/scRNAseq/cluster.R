set.seed(1234) # optional
library(Seurat)
library(clustree)
library(reticulate)

sample <- "ATsim6d"
# 读取
(obj <- readRDS(paste0(file.path("./Seurat_obj", sample), ".rds")))

obj <- NormalizeData(
  object = obj,
  normalization.method = "LogNormalize",
  scale.factor = 1e4
)

# selection.method 有三种方法：vst，mean.var.plot，dispersion。
# 默认选择2000个HVG……咱也不知道哪种方法合适
# TQ师兄、陈瑜用mean.var.plot
obj <- FindVariableFeatures(
  object = obj,
  selection.method = "mean.var.plot",
  mean.cutoff = c(0.0125, 3),
  dispersion.cutoff = c(1.0, Inf)
)
# 细胞周期score
cc_gene <- read.csv("/data5/wmc_data/reference/TAIR10/Arabidopsis_cell_cycle_genes_fx.csv")
g1s <- cc_gene[cc_gene$phase == "G1S", 1]
g2m <- cc_gene[cc_gene$phase == "G2M", 1]

obj <- CellCycleScoring(
  obj,
  s.features = g1s,
  g2m.features = g2m,
  set.ident = TRUE
)

obj[["CC.Difference"]] <- obj$S.Score - obj$G2M.Score

# 很慢
obj <- ScaleData(
  object = obj,
  features = rownames(x = obj),
  vars.to.regress = c("nCount_RNA", "percent.mito", "percent.chlo", "CC.Difference")
)

## 线性降维 -------
obj <- RunPCA(
  object = obj,
  features = VariableFeatures(obj),
  verbose = FALSE,
  npcs = 100
)

# 判断dim选多少
ElbowPlot(obj, reduction = "pca", ndims = 100)

dimension <- c(1:60)
res <- seq(.1, 1.6, .2)

obj <- FindNeighbors(obj, dims = dimension)
# resolution 可以多次调试的,不同resolution的分群
tmp <- FindClusters(
  object = obj,
  resolution = res
)

clustree(tmp@meta.data, prefix = "RNA_snn_res.")
rm(tmp)
gc()

obj <- FindClusters(
  object = obj,
  resolution = 0.8
)

obj <- RunTSNE(obj, dims = dimension)
obj <- RunUMAP(
  object = obj,
  reduction = "pca",
  umap.method = "uwot",
  dims = dimension,
  n.neighbors = 30,
  metric = "correlation",
  min.dist = 0.3
)


DimPlot(obj, reduction = "tsne", label = TRUE, pt.size = 1)

p_umap <- DimPlot(obj, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1) +
  tidydr::theme_dr() +
  theme(panel.grid = element_blank())

ggsave(file.path("Plots", paste0(sample, "_umap.png")), p_umap, width = 10, height = 8)

saveRDS(obj, file.path("./Seurat_obj", paste0(sample, ".rds")))
