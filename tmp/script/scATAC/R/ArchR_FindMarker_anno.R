library(parallel)
library(ArchR)
addArchRThreads(threads = 32)

proj <- loadArchRProject("Result/BasalStem-save")
# identify marker genes with genescore
markersGS <- getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "GeneScoreMatrix",
    groupBy = "Clusters",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1")

# 输出结果 -----
alias <- read.csv("~/reference/Oryza_sativa/annotation/IRGSP_annotation_At_simple.csv")
alias <- alias[!duplicated(alias$RAP_ID), ]

markerList <- lapply(markerList, function(x) {
    if (nrow(x) < 1) {
        return(x)
    }
    x$Description <- ""
    x$Symbol <- ""
    x$TAIR <- ""
    for (i in seq_len(nrow(x))) {
        gene <- x$name[i]
        if (gene %in% alias$RAP_ID) {
            des <- alias[alias$RAP_ID %in% gene, 2]
            symbol <- alias[alias$RAP_ID %in% gene, 3]
            tair <- alias[alias$RAP_ID %in% gene, 4]
            x$Description[i] <- des
            x$Symbol[i] <- symbol
            x$TAIR[i] <- tair
        }
    }
    x
})

openxlsx::write.xlsx(markerList, "MarkerGene/BasalStem_Cluster_MarkerGene.xlsx")


p <- plotEmbedding(
    ArchRProj = proj,
    colorBy = "GeneScoreMatrix",
    name = "Os02g0537775",
    pal = ArchRPalettes$whiteBlue,
    embedding = "UMAP",
    rastr = FALSE,
    plotAs = "point",
    quantCut = c(0.01, 0.99),
    imputeWeights = NULL
)
p


# GO enrichment --------

library(clusterProfiler)
library(org.At.tair.db)

ego_bp <- enrichGO(markerList$C1$name, org.At.tair.db, keyType = "TAIR", ont = "MF")
dotplot(ego_bp)

# plot----
plotEmbedding(
    ArchRProj = proj,
    colorBy = "cellColData",
    name = "Clusters",
    embedding = "UMAP"
)


df <- getCellColData(proj, select = c("log10(nFrags)", "TSSEnrichment"))

ggPoint(
    x = df[, 1],
    y = df[, 2],
    colorDensity = TRUE,
    continuousSet = "sambaNight",
    xlabel = "Log10 Unique Fragments",
    ylabel = "TSS Enrichment",
    xlim = c(log10(3000), quantile(df[, 1], probs = 0.99)),
    ylim = c(1.2, quantile(df[, 2], probs = 0.99))
) + geom_hline(yintercept = 1.5, lty = "dashed") +
    geom_vline(xintercept = 3.5, lty = "dashed")

plotGroups(
    ArchRProj = proj,
    groupBy = "Sample",
    colorBy = "cellColData",
    name = "TSSEnrichment",
    plotAs = "ridges"
)

plotFragmentSizes(ArchRProj = proj)

plotTSSEnrichment(ArchRProj = proj)

table(proj$Clusters)


heatmapGS <- plotMarkerHeatmap(
    seMarker = markersGS,
    cutOff = "FDR <= 0.01 & Log2FC >= 1.25",
    transpose = TRUE
)
# pdf("tmp.pdf", width = 40, height = 12)
# ComplexHeatmap::draw(
#     heatmapGS,
#     heatmap_legend_side = "bot"
#     # annotation_legend_side = "bot"
# )
# dev.off()
plotPDF(heatmapGS,
    name = "GeneScores-Marker-Heatmap",
    width = 40, height = 12,
    ArchRProj = proj, addDOC = FALSE
)


cM <- confusionMatrix(paste0(proj$cellNames), paste0(proj$Clusters))

cM <- cM / Matrix::rowSums(cM)