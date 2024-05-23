library(Seurat)
h5_files <- list.files("data/", full.names = TRUE) 
# load DS data
seurat.objs <- list()
for (sample in h5_files) {
  sample_ids <- strsplit(basename(sample), "_", fixed=TRUE)[[1]]
  sample_id <- sample_ids[1]
  mtx <-  Read10X_h5(sample)
  seurat.objs[[sample_id]] <- CreateSeuratObject(counts = mtx, 
                                                 project = sample_id, 
                                                 min.cells = 3, 
                                                 min.features = 200)
  
}
#merge together, combine
seurat.obj <- merge(x=seurat.objs[[1]], y = seurat.objs[-1])
rm(seurat.objs);gc()
#---------------------------------------------------------
seurat.obj[["percent.ch"]] <- PercentageFeatureSet(seurat.obj,pattern = "^ATCG")
seurat.obj[['ratio']] <- seurat.obj$nCount_RNA / seurat.obj$nFeature_RNA
m <- mean(seurat.obj$ratio) 
s <- 2 * sd(seurat.obj$ratio)

seurat.obj[['outlier']] <- ifelse(seurat.obj$ratio > m + s |  
                                    seurat.obj$ratio < m - s, "out", "in")

# filtration
seurat.obj <- subset(seurat.obj, subset = percent.ch < 40 & outlier == "in")
seurat.obj

# Normalization ------------------------------------
seurat.obj <- NormalizeData(seurat.obj, normalization.method = "LogNormalize", 
                          scale.factor = 10000)
seurat.obj <- FindVariableFeatures(seurat.obj, selection.method = "vst", nfeatures = 2000)
all.genes<- rownames(seurat.obj)
seurat.obj <- ScaleData(seurat.obj, features = all.genes)

# linear dimensional reduction -----------------------------------------------
seurat.obj <- RunPCA(seurat.obj, features = VariableFeatures(object = seurat.obj),npcs = 100)
DimPlot(seurat.obj, reduction = "pca")
ElbowPlot(seurat.obj,reduction = "pca",ndims = 50)

seurat.obj <- FindNeighbors(seurat.obj, reduction = "pca", dims = 1:50)
seurat.obj <- FindClusters(seurat.obj, resolution = 0.21)
seurat.obj <- RunUMAP(seurat.obj, reduction = "pca", dims = 1:50,
                      n.neighbors = 30,
                      min.dist = 0.3)
#seurat.obj <- RunTSNE(seurat.obj, reduction = "harmony", dims = 1:50)

DimPlot(seurat.obj, reduction = "umap",split.by = "orig.ident")
DimPlot(seurat.obj,reduction = "umap")
FeaturePlot(seurat.obj,features = c("AT3G26744"),split.by = "orig.ident")


