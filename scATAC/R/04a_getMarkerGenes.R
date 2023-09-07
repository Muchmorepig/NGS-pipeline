library(parallel)
suppressMessages(library(ArchR))
library(fs)
addArchRThreads(threads = 32)

# load ArchR project
proj_path <- "Result/SIM12d"
proj <- loadArchRProject(proj_path)

# identify marker genes with genescore -------
markersGS <- getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "GeneScoreMatrix",
    groupBy = "Tile_4",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)


heatmapGS <- plotMarkerHeatmap(
    seMarker = markersGS,
    cutOff = "FDR <= 0.01 & Log2FC >= 1",
    transpose = TRUE
)

plotdir <- paste0(proj_path, "/Plots")
if (!dir_exists(plotdir)) dir_create(plotdir)

pdf(paste0(plotdir, "/markerHeatmap.pdf"), width = 60, height = 12)
ComplexHeatmap::draw(
    heatmapGS,
    heatmap_legend_side = "bot"
    # annotation_legend_side = "bot"
)
dev.off()


markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1")

# 输出结果 -----
# library(data.table)
alias <- fread("~/reference/Oryza_sativa/annotation/IRGSP_OStair_anno.csv")
table(duplicated(alias))

markerPeaksList <- lapply(markerPeaksList, function(x) {
  if (nrow(x) < 1) {
    return(x)
  }
  sub_anno <- alias[RAP_ID %in% x$RAP_ID, ]
  x <- dplyr::left_join(x, sub_anno)
  return(x)
})

markerdir <- paste0(proj_path, "/Marker")
if (!dir_exists(markerdir)) dir_create(markerdir)

openxlsx::write.xlsx(markerList, paste0(markerdir, "/MarkerGene_Tile4.xlsx"))


#####################
# 测试-------
library(magrittr)
iter2 <- markerList
iter2 <- sapply(seq.int(iter2), function(x) {
    iter2[[x]]$Cluster <- names(iter2[x])
    iter2[[x]]
}) %>% purrr::reduce(., rbind)

library(ggplot2)
dat <- data.frame(
    iterTimes = c("2", "4", "8"),
    Markers_Number = c(nrow(iter2), nrow(iter4), nrow(iter8))
)

df <- data.frame(
    "per" = rep(c("dup", "uniq"), 3),
    "iterTimes" = c("2", "2", "4", "4", "8", "8"),
    "Marker_Number" = c(170, 2656, 578, 3381, 879, 3886)
)

ggplot(df, aes(x = iterTimes, y = Marker_Number, fill = per)) +
    geom_bar(stat = "identity") +
    scale_fill_brewer(palette = "Paired") +
    theme_minimal()


a <- unique(iter2$name)
b <- unique(iter4$name)
d <- unique(iter8$name)
e <- intersect(a, b)
f <- intersect(b, d)
g <- intersect(a, d)

h <- intersect(e, d)
i <- length(e) - length(h)
j <- length(g) - length(h)
k <- length(f) - length(h)

length(a) - i - length(h) - j
length(b) - i - length(h) - k
length(d) - j - length(h) - k

library(d3vennR)
vn <- d3vennR(
    data = list(
        list(sets = list(0), label = "A 506", size = length(a) / 2000),
        list(sets = list(1), label = "B 840", size = length(b) / 2000),
        list(sets = list(2), label = "D 1225", size = length(d) / 2000),
        list(sets = c(0, 1), label = "I 207", size = length(e) / 2000),
        list(sets = c(0, 2), label = "J 327", size = length(g) / 2000),
        list(sets = c(1, 2), label = "K 718", size = length(f) / 2000),
        list(sets = c(0, 1, 2), label = "H 1616", size = length(h) / 2000)
    )
)

htmlwidgets::saveWidget(vn,
    file = "Result/BasalStem/Plots/vn.html",
    selfcontained = FALSE
)