suppressMessages(library(ArchR))
library(parallel)
library(magrittr)
library(fs)
addArchRThreads(threads = 32)
addArchRChrPrefix(chrPrefix = FALSE)
load("MSU7.Rdata")

set.seed(1234) # 可加可不加

arrowFiles <- c(
    # "Result/BasalStem/ArrowFiles/BasalStem.arrow",
    "Result/CIM22d/ArrowFiles/CIM22d.arrow",
    "Result/SIM12d/ArrowFiles/SIM12d.arrow",
    "Result/SIM16d/ArrowFiles/SIM16d.arrow"
)

output_dir <- "Result/Merge"
if (!dir_exists(output_dir)) dir_create(output_dir)


proj <- ArchRProject(
    ArrowFiles = arrowFiles,
    outputDirectory = output_dir,
    copyArrows = TRUE,
    geneAnnotation = geneAnnotation,
    genomeAnnotation = genomeAnnotation,
)

iter <- 4
lsi_name <- "iter4"

proj %<>%
    addIterativeLSI(
        ArchRProj = .,
        useMatrix = "TileMatrix",
        name = lsi_name,
        iterations = iter,
        clusterParams = list(
            resolution = c(0.8),
            sampleCells = 10000,
            n.start = 10
        ),
        force = TRUE
    ) %>%
    addClusters(
        input = .,
        reducedDims = lsi_name,
        name = paste0("Cluster_", lsi_name),
        force = TRUE
    ) %>%
    addUMAP(
        ArchRProj = .,
        reducedDims = lsi_name,
        nNeighbors = 40,
        minDist = 0.4,
        metric = "cosine",
        name = paste0("UMAP_", lsi_name),
        force = TRUE
    )


# 保存项目
saveArchRProject(proj)
