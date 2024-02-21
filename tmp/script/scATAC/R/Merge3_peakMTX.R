suppressMessages(library(ArchR))
library(parallel)
library(magrittr)
library(fs)
addArchRThreads(threads = 32)

proj_path <- "Result/Merge"
proj <- loadArchRProject(proj_path)

proj <- addGroupCoverages(
    ArchRProj = proj,
    groupBy = "Cluster_iter4",
    minCells = 15,
    maxCells = 1500,
    force = TRUE
)

## call peaks
proj <- addReproduciblePeakSet(
    ArchRProj = proj,
    groupBy = "Cluster_iter4",
    genomeSize = 3.7e08,
    excludeChr = c("ChrM", "ChrC"),
    pathToMacs2 = "/opt/biosoft/MACS2.1.2/bin/macs2",
    force = TRUE
)

## calculate peak matrix
proj <- addPeakMatrix(proj, binarize = TRUE, force = TRUE)

proj %<>% addIterativeLSI(
    ArchRProj = .,
    useMatrix = "PeakMatrix",
    name = "Peaks",
    varFeatures = 25000,
    iterations = 4,
    force = TRUE
) %>%
    addClusters(
        input = .,
        reducedDims = "Peaks",
        name = "Cluster_Peaks",
        resolution = 0.8,
        force = TRUE,
    ) %>%
    addUMAP(
        ArchRProj = .,
        reducedDims = "Peaks",
        name = "UMAP_Peaks",
        nNeighbors = 30,
        force = TRUE
    )

saveArchRProject(proj)

source("./bin/umapForArchR.R")

umap_clu <- umapForArchR(proj, embed.umap = "UMAP_Peaks",
                     embed.cluster = "Cluster_Peaks",
                     point.size = 1.2,
                     theme = "clear",
                     title = "UMAP Clustered by PeakMatrix")
umap_sample <- umapForArchR(proj, embed.umap = "UMAP_Peaks",
                         embed.cluster = "Sample",
                         point.size = 1.2,
                         theme = "clear",
                         title = "UMAP Clustered by PeakMatrix")

dir_exists(paste0(proj_path, "/Plots"))

library(patchwork)

pdf(paste0(proj_path, "/Plots/merge_umap_peaks.pdf"), width = 14)
umap_sample + umap_clu
dev.off()