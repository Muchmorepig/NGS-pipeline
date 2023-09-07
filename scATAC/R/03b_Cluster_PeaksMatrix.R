suppressMessages(library(ArchR))
library(parallel)
library(magrittr)
addArchRThreads(threads = 32)

set.seed(1234) # 可加可不加
# load ArchR project
proj_path <- "Result/SIM16d" # 之前的 output dir
proj <- loadArchRProject(path = proj_path)

## create pseudo-bulk replicates
proj <- addGroupCoverages(
  ArchRProj = proj,
  groupBy = "Tile_4",
  minCells = 15,
  maxCells = 1200,
  force = TRUE
)

## call peaks
proj <- addReproduciblePeakSet(
  ArchRProj = proj,
  groupBy = "Tile_4",
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
) %>% addClusters(
  input = .,
  reducedDims = "Peaks",
  name = "Cluster_Peaks",
  resolution = 0.8
) %>% addUMAP(
  ArchRProj = .,
  reducedDims = "Peaks",
  name = "UMAP_Peaks",
  nNeighbors = 40,
  force = TRUE
)

saveArchRProject(proj)

source("./bin/umapForArchR.R")

umap <- umapForArchR(proj, embed.umap = "UMAP_Peaks", embed.cluster = "Cluster_Peaks", point.size = 1.8, title = "UMAP of SIM16d (Clustered by PeakMatrix)")
umap2 <- umapForArchR(proj, embed.umap = "Tile_4", embed.cluster = "Tile_4", point.size = 1.8, title = "UMAP of SIM16d (Clustered by TileMatrix)")

library(patchwork)

pdf(paste0(proj_path, "/Plots/umap.pdf"), width = 14)
umap2 + umap
dev.off()

