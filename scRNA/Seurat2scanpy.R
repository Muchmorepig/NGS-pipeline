library(Seurat)
library(SeuratDisk)

obj <- readRDS("./scRNA/Seurat_obj/Merge_strata.rds")

SaveH5Seurat(obj, filename = "Merge_strata.h5Seurat")
Convert("Merge_strata.h5Seurat", dest = "h5ad")
