library(parallel)
suppressMessages(library(ArchR))
addArchRThreads(threads = 32)

# load ArchR project
proj_path <- "Result/CIM22d" # 之前的 output dir
proj <- loadArchRProject(path = proj_path)

# load scRNA-seq seurat object
seRNA <- readRDS("rawdata/scRNA_CIM22d_CY.rds")

getAvailableMatrices(proj)

proj <- addGeneIntegrationMatrix(
    ArchRProj = proj,
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    # reducedDims = "Tile_4",
    reducedDims = "Peaks",
    seRNA = seRNA,
    addToArrow = FALSE,
    # groupATAC = "Tile_4",
    groupATAC = "Cluster_Peaks",
    groupRNA = "seurat_clusters",
    # nameCell = "predictedCell_Un",
    nameCell = "predictedCell_Un2",
    # nameGroup = "predictedGroup_Un",
    nameGroup = "predictedGroup_Un2",
    # nameScore = "predictedScore_Un",
    nameScore = "predictedScore_Un2",
    force = TRUE
)


cM2 <- as.matrix(confusionMatrix(
    proj$Cluster_Peaks,
    proj$predictedGroup_Un2
))

col_fun <- circlize::colorRamp2(c(0, max(log10(cM + 1))), c("white", "#f54434"))
ComplexHeatmap::Heatmap(
    column_title = "Log10(Cell Number + 1)",
    log10(cM2 + 1),
    show_heatmap_legend = FALSE,
    col = col_fun,
    cluster_rows = FALSE,
    # cluster_columns = F,
    row_names_side = "left"
)

preClust <- colnames(cM)[apply(cM, 1, which.max)]
preClust2 <- colnames(cM2)[apply(cM2, 1, which.max)]

cbind(preClust, rownames(cM))
cbind(preClust2, rownames(cM2)) # Assignments

# peakmat <- getMatrixFromProject(proj, "PeakMatrix", binarize = TRUE)
