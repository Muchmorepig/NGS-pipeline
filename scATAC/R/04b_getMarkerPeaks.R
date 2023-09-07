suppressMessages(library(ArchR))
library(parallel)
library(magrittr)
addArchRThreads(threads = 32)

proj_path <- "./Result/SIM12d"
proj <- loadArchRProject(proj_path)
getAvailableMatrices(proj)

## 获取marker Peaks -------
markerPeaks <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "PeakMatrix",
  groupBy = "Cluster_Peaks",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

(markerPeaksList <- getMarkers(
  markerPeaks,
  cutOff = "FDR <= 0.05 & Log2FC >= 1"
))
sapply(markerPeaksList, nrow)

# Peaks注释 -------------------
# library(FindIT2)
library(TxDb.Osativa.MSU.xzg)
txdb <- TxDb.Osativa.MSU.xzg
seqlevels(txdb)

markerPeaksList <- lapply(markerPeaksList, function(x) {
  # print(x)
  if (nrow(x) < 1) {
    return(x)
  }

  peak_gr <- GRanges(
    seqnames = x$seqnames,
    IRanges(start = x$start, end = x$end),
    feature_id = paste0("peak_", row.names(x))
  )
  mm_anno_ng <- FindIT2::mm_nearestGene(
    peak_GR = peak_gr,
    Txdb = txdb
  )
  x$RAP_ID <- mm_anno_ng$gene_id
  x$distanceToTSS <- mm_anno_ng$distanceToTSS
  x
})

# 输出结果
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
if (!dir.exists(markerdir)) dir.create(markerdir)


openxlsx::write.xlsx(markerPeaksList, paste0(markerdir, "/MarkerPeaks_iter4.xlsx"))
