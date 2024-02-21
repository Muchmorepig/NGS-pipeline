library(Seurat)
library(magrittr)
# 创建Seurat对象
sample <- c(
  ATsim4d = "./ATSIM4D/outs/filtered_feature_bc_matrix"
)

(scrna <- Read10X(data.dir = sample) %>%
  CreateSeuratObject(.,
    min.cells = 3,
    min.features = 200,
    project = names(sample)
  ))

# 添加细胞器基因比例
scrna[["percent.organelle"]] <- PercentageFeatureSet(scrna, pattern = "^ATCG|^ATMG")
scrna[["percent.mito"]] <- PercentageFeatureSet(scrna, pattern = "^ATMG")
scrna[["percent.chlo"]] <- PercentageFeatureSet(scrna, pattern = "^ATCG")

# 看一下，考虑一下过滤条件
source("/data5/wmc_data/utils/qc_plot.R")

(p <- qc_plot(scrna, max_nFeature = 8000, min_nFeature = 750))

if (!dir.exists("./Plots")) dir.create("./Plots", recursive = TRUE)
ggsave(paste0("./Plots/", names(sample), "_qc.pdf"), p, width = 10, height = 8)

# 考虑低表达基因
# sel_genes <- rownames(obj)[Matrix::rowSums(obj@assays$RNA@counts>0) > 3]

# 过滤数据
(scrna <- subset(scrna,
  subset = nFeature_RNA > 750 &
    nFeature_RNA < 8000 &
    percent.mito < 5 & percent.chlo < 5
  # features = sel_genes
))

# 细胞很多的话可以去除双细胞
# 保存
saveRDS(scrna, paste0("./Seurat_obj/", names(sample), ".rds"))
