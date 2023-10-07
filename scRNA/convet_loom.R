library(Seurat)
library(SeuratDisk)

obj <- readRDS("./integrated/callusSAM_stage4_9samples/seurat4/code/callusSAM_S9_S4AtR27_singlet_harmony_beta1.Rds")

obj <- subset(obj, subset = seurat_clusters %in% c(0 ,18, 19))
table(obj$stage)

obj.loom <- as.loom(obj, filename = "./loom/c_0_18_19.loom", verbose = FALSE)
obj.loom
obj.loom$close_all()