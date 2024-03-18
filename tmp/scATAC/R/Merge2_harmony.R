suppressMessages(library(ArchR))
library(parallel)
library(magrittr)
library(fs)
addArchRThreads(threads = 32)

proj_path <- "Result/Merge"
proj <- loadArchRProject(proj_path)

proj <- addHarmony(
    ArchRProj = proj,
    reducedDims = "iter4",
    name = "HarmonyIter4",
    groupBy = "Sample",
    force = TRUE
)

proj %<>% addClusters(
    input = .,
    reducedDims = "HarmonyIter4",
    name = "ClusterHarmony_iter4",
    force = TRUE
) %>% addUMAP(
    ArchRProj = .,
    reducedDims = "HarmonyIter4",
    name = "UMAP_Harmony_iter4",
    nNeighbors = 40,
    minDist = 0.4,
    metric = "cosine",
    force = TRUE
)

saveArchRProject(proj)

source("./bin/umapForArchR.R")

p1 <- umapForArchR(proj,
    embed.cluster = "Sample",
    embed.umap = "UMAP_iter4",
    title = "Colored By Sample",
    point.size = 1.2
)

p2 <- umapForArchR(proj,
    embed.cluster = "Cluster_iter4",
    embed.umap = "UMAP_iter4",
    title = "Colored By Cluster",
    point.size = 1.2
)

p3 <- umapForArchR(proj,
    embed.cluster = "Sample",
    embed.umap = "UMAP_Harmony_iter4",
    title = "Harmony Colored By Sample",
    point.size = 1.2
)

p4 <- umapForArchR(proj,
    embed.cluster = "ClusterHarmony_iter4",
    embed.umap = "UMAP_Harmony_iter4",
    title = "Harmony Colored By Cluster",
    point.size = 1.2
)

library(patchwork)


pdf("./Result/Merge/Plots/merge_umap.pdf", height = 12, width = 12)
(p1 + p2) / (p3 + p4)
dev.off()
