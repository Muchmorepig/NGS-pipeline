library(ArchR)
addArchRThreads(threads = 32)
# 加载项目 ------
proj <- "Result/Merge"
proj <- loadArchRProject(proj)
# 查看项目信息 -------
colnames(proj@cellColData)
getAvailableMatrices(proj)

# 获取 Marker Peaks -------
markerPeaks <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "PeakMatrix",
  groupBy = "Cluster_Peaks",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
# 设置阈值过滤 FDR <= 0.05 & Log2FC >= 1
markerPeaksList <- getMarkers(
  markerPeaks,
  cutOff = "FDR <= 0.05 & Log2FC >= 1"
)

sapply(markerPeaksList, nrow)

# Marker Peaks注释 ------
library(TxDb.Osativa.MSU.xzg)
txdb <- TxDb.Osativa.MSU.xzg
seqlevels(txdb)

# Peak -> gene_id & distanceToTSS
markerPeaksList <- lapply(markerPeaksList, function(x) {
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
  as.data.frame(x)
})

# 添加 symbol、desscription、TAIR 信息
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

# 保存结果 ------
markerdir <- paste0(getOutputDirectory(proj), "/Marker")
if (!dir.exists(markerdir)) dir.create(markerdir)

openxlsx::write.xlsx(markerPeaksList, paste0(markerdir, "/MarkerPeaks.xlsx"))
