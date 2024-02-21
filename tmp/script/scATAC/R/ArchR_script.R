# 加载R包和函数
source("scripts/utilize.R")
library(openxlsx)
library(ArchR)

addArchRThreads(threads = 32)
addArchRChrPrefix(chrPrefix = FALSE)
load("data/ref/TAIR10.Rdata")


# 预处理 ---------------------------------------------------------------------
frag_file <- "./data/fragments/sim4d_fragments.tsv.gz"
sample_name <- "sim4d_raw"
arrowFiles <- createArrowFiles(
  inputFiles = frag_file,
  sampleNames = sample_name,
  geneAnnotation = geneAnnotation,
  genomeAnnotation = genomeAnnotation,
  minTSS = 0,
  minFrags = 100,
  maxFrags = 1e6,
  promoterRegion = c(1000, 500),
  excludeChr = c("ChrM", "ChrC"),
  addTileMat = FALSE,
  addGeneScoreMat = FALSE,
  force = TRUE
)

## 将在当前目录下生成arrow移动到给定目录下
if (!dir.exists("data/raw_arrows")) {
  dir.create("data/raw_arrows", recursive = TRUE)
}

arrowFiles <- move_file(arrowFiles, "data/raw_arrows")


output_dir <- "temp" # 输出文件夹
proj <- ArchRProject(
  ArrowFiles = arrowFiles,
  outputDirectory = output_dir,
  copyArrows = TRUE, # This is recommened so that you maintain an unaltered copy for later usage.,
  geneAnnotation = geneAnnotation,
  genomeAnnotation = genomeAnnotation,
)

df <- getCellColData(proj, select = c("log10(nFrags)", "TSSEnrichment", "Sample"))
p <- ggPoint(
  x = df[, 1],
  y = df[, 2],
  colorDensity = TRUE,
  continuousSet = "sambaNight",
  size = 0.75,
  xlabel = "Log10 Unique Fragments(SIM 0d)",
  ylabel = "TSS Enrichment",
  xlim = c(log10(100), log10(100000)),
  ylim = c(0, 4),
  rastr = TRUE
) + geom_hline(yintercept = 1.5, lty = "dashed") + geom_vline(xintercept = log10(1000), lty = "dashed")


cell_number <- nrow(df[df$`log10(nFrags)` > log10(1000) & df$TSSEnrichment > 1.5, ])

p <- p + ggtitle(paste0("SIM 4D: cell number ", cell_number))

# 保存输出结果

ggsave("results/figures/sim4d_qc.pdf", p)


# 正式分析 --------------------------------------------------------------------

## 设置过滤参数
minTSS <- 1.5
minFrags <- 1000
maxFrags <- 1e5

frag_file <- "./data/fragments/sim4d_fragments.tsv.gz"
samle_name <- "sim4d"
## 创建arrow文件，包含gene score matrix
arrowFiles <- createArrowFiles(
  inputFiles = frag_file,
  sampleNames = samle_name,
  geneAnnotation = geneAnnotation,
  genomeAnnotation = genomeAnnotation,
  minTSS = minTSS,
  minFrags = minFrags,
  maxFrags = maxFrags,
  promoterRegion = c(2000, 1000),
  excludeChr = c("ChrM", "ChrC"),
  addTileMat = TRUE,
  TileMatParams = list(tileSize = 100),
  addGeneScoreMat = TRUE,
  force = TRUE
)
## 计算doublet score, 用于后续过滤潜在的doublet
doubScores <- addDoubletScores(
  input = arrowFiles,
  k = 10,
  knnMethod = "UMAP",
  LSIMethod = 1
)

## 将在当前目录下生成arrow移动到给定目录下
if (!dir.exists("data/raw_arrows")) {
  dir.create("data/raw_arrows", recursive = TRUE)
}

arrowFiles <- file.move(arrowFiles, "data/raw_arrows")


##  设置项目输出文件夹
output_dir <- "results/SIM4D"

# 创建 ArchRProject
proj <- ArchRProject(
  ArrowFiles = arrowFiles,
  outputDirectory = output_dir,
  copyArrows = TRUE, # This is recommened so that you maintain an unaltered copy for later usage.,
  geneAnnotation = geneAnnotation,
  genomeAnnotation = genomeAnnotation,
)
## 过滤doublet
proj <- filterDoublets(proj)

## 使用TileMatrix以迭代LSI的方式进行降维
### 主要修改iterations和clusterParams
proj <- addIterativeLSI(
  ArchRProj = proj,
  useMatrix = "TileMatrix",
  name = "IterativeLSI",
  iterations = 3,
  clusterParams = list(
    resolution = c(2, 1),
    sampleCells = 10000,
    maxClusters = c(6, 10)
  ),
  dimsToUse = 1:50,
  force = TRUE
)

## 聚类分析
### 默认方法是Seurat, 设置resolution可以修改类群数目
proj <- addClusters(
  input = proj,
  reducedDims = "IterativeLSI",
  name = "Clusters",
  resolution = 0.8,
  force = TRUE
)

## 可以按照Cluster做一些基本统计
library(dplyr)
as_tibble(proj@cellColData) %>%
  group_by(Clusters) %>%
  summarise(min(TSSEnrichment), mean(TSSEnrichment), median(TSSEnrichment), max(TSSEnrichment)) %>%
  ungroup()


## 添加UMAP非线性降维
proj <- addUMAP(
  ArchRProj = proj,
  reducedDims = "IterativeLSI",
  nNeighbors = 40,
  minDist = 0.4,
  metric = "cosine",
  force = TRUE
)

## 保存ArchRProject的元信息
## 这一行需要随时调用，确保数据不丢失
saveRDS(proj, file = file.path(output_dir, "ArchRProject.Rds"))

## 加载数据
proj <- readRDS(file = file.path(output_dir, "ArchRProject.Rds"))

## UMAP可视化部分

## plotEmbedding: colorBy 和 name 用于获取分组

## 展示Cluster
## cellColdata是记录分组的元信息, Clusters是其中的一个列名
p <- plotEmbedding(
  ArchRProj = proj,
  colorBy = "cellColData",
  name = "Clusters", embedding = "UMAP"
)
ggsave(file.path(output_dir, "UMAP_by_Cluaters.pdf"), p)

## 展示readsInTSS
p <- plotEmbedding(
  ArchRProj = proj,
  colorBy = "cellColData",
  name = "ReadsInTSS", embedding = "UMAP"
)
ggsave(file.path(output_dir, "UMAP_by_ReadsInTSS.pdf"), p)

p <- plotEmbedding(
  ArchRProj = proj,
  colorBy = "cellColData",
  name = "TSSEnrichment", embedding = "UMAP"
)
ggsave(file.path(output_dir, "UMAP_by_TSSEnrichment.pdf"), p)


# 差异分析 --------------------------------------------------------------------


## 基于genescore ----------------------------------------------------------
markersGS <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "GeneScoreMatrix",
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

## FDR(adjust p.value) Log2FC(log2 Fold Change)
markersGSList <- getMarkers(markersGS, cutOff = "FDR <= 0.05 & Log2FC >= 1")

## 导出差异基因
### 增加别名
alias <- read.table("/data/reference/genome/TAIR10/ath_gene_alias.txt", fill = TRUE, stringsAsFactors = FALSE)
alias <- alias[!duplicated(alias$V1), ]

markersGSList <- lapply(markersGSList, function(x) {

  # print(x)
  if (nrow(x) < 1) {
    return(x)
  }
  x$alias <- ""
  for (i in seq.int(1, nrow(x))) {
    gene <- x$name[i]
    if (gene %in% alias$V1) {
      gene_symbol <- alias[alias$V1 %in% gene, 2]
      x$alias[i] <- gene_symbol
    }
  }
  x
})

## 导出为xlsx
write.xlsx(markersGSList, file.path(output_dir, "markersGSList_FDR0.05_LogFC1.xlsx"))

## 差异基因展示
## pal是调色盘
p <- plotEmbedding(
  ArchRProj = proj,
  colorBy = "GeneScoreMatrix",
  name = "AT1G75240",
  embedding = "UMAP",
  quantCut = c(0.01, 0.95),
  imputeWeights = NULL,
  rastr = FALSE,
  pal = ArchRPalettes$whiteBlue,
  # continuousSet = p
  plotAs = "points"
)
p
## TODO List:
### 增加小提琴图，增加UMAP的cluster标签

## 输出差异基因的热图
heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS,
  cutOff = "FDR <= 0.05 & Log2FC >= 1",
  # labelMarkers = markerGenes,
  nLabel = 3,
  transpose = FALSE
)
heatmapGS

## 基于peak   ----------------------------------------------------------

## 先找到每个cluster的peak
###  MACS2的路径
pathToMacs2 <- "/opt/biosoft/MACS2.1.2/bin/macs2"

### 添加 group coverage 用于标准化
### proj@projectMetadata$GroupCoverages
proj <- addGroupCoverages(proj)

### 运行MACS2 按照clusters寻找peak
proj <- addReproduciblePeakSet(
  ArchRProj = proj,
  groupBy = "Clusters",
  pathToMacs2 = pathToMacs2,
  genomeSize = "110000000",
  excludeChr = c("ChrM", "ChrC"),
  genomeAnnotation = genomeAnnotation,
  geneAnnotation = geneAnnotation
)

## 将peak加入到arrow中, 构建peak x cell的count矩阵
## 保存为 PeakMatrix, 可用于cisTopic的分析
proj <- addPeakMatrix(proj)

saveRDS(proj, file = file.path(output_dir, "ArchRProject.Rds"))

## 查看proj当前的matrix
getAvailableMatrices(proj)

## 差异peak开放状态分析
markersPeaks <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "PeakMatrix",
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markersPeaksList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.05 & Log2FC >= 1")
markersPeaksList

## 对peak进行注释
library(FindIT2)
library(TxDb.Athaliana.BioMart.plantsmart28)
txdb <- TxDb.Athaliana.BioMart.plantsmart28
seqlevels(txdb) <- paste0("Chr", c(1:5, "M", "C"))

markersPeaksList <- lapply(markersPeaksList, function(x) {

  # print(x)
  if (nrow(x) < 1) {
    return(x)
  }

  peak_gr <- GRanges(seqnames = x$seqnames, IRanges(
    start = x$start, end = x$end
  ), feature_id = paste0("peak_", row.names(x)))
  mm_anno_ng <- mm_nearestGene(
    peak_GR = peak_gr,
    Txdb = txdb
  )
  x$gene_id <- mm_anno_ng$gene_id
  x$distanceToTSS <- mm_anno_ng$distanceToTSS
  x
})

## 为了方便阅读，添加别名
alias <- read.table("/data/reference/genome/TAIR10/ath_gene_alias.txt", fill = TRUE, stringsAsFactors = FALSE)
alias <- alias[!duplicated(alias$V1), ]

markersPeaksList <- lapply(markersPeaksList, function(x) {

  # print(x)
  if (nrow(x) < 1) {
    return(x)
  }
  x$alias <- ""
  for (i in seq.int(1, nrow(x))) {
    gene <- x$gene_id[i]
    if (gene %in% alias$V1) {
      gene_symbol <- alias[alias$V1 %in% gene, 2]
      x$alias[i] <- gene_symbol
    }
  }
  x
})

## 导出导出差异peak的xlsx
write.xlsx(markersGSList, file.path(output_dir, "markersGSList_FDR0.05_LogFC1.xlsx"))


## 导出peak 和 bigwig 用于IGV
peakset <- getPeakSet(proj)

peak_df <- data.frame(
  seqnames = seqnames(peakset),
  start = start(peakset), end = end(peakset),
  score = peakset$score, strand = "*", cluster = names(peakset)
)
peak_bed_path <- file.path(output_dir, "scATAC_peak.bed")
rtracklayer::export.bed(peak_df, peak_bed_path)

getGroupBW(proj, groupBy = "Clusters", normMethod = "ReadsInTSS")


# Motif富集分析 ---------------------------------------------------------------

## 添加Motif信息
proj <- addMotifAnnotations(
  ArchRProj = proj,
  motifSet = "JASPAR2020", name = "Motif",
  species = "Arabidopsis thaliana"
)
## 分析每个cluster的peak的富集motif
enrichMotifs <- peakAnnoEnrichment(
  seMarker = markersPeaks,
  ArchRProj = proj,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)

# 伪时间分析 -------------------------------------------------------------------

# trajectory <- paste0("C", c(1:2, 5:8) )
trajectory <- "C7"
proj <- addTrajectory(
  ArchRProj = proj,
  name = "C7",
  groupBy = "Clusters",
  trajectory = trajectory,
  embedding = "UMAP",
  force = TRUE
)


p <- plotTrajectory(proj,
  trajectory = "C7",
  colorBy = "cellColData",
  name = "C7"
)
p[[1]]
p <- plotTrajectory(projHeme5,
  trajectory = "MyeloidU",
  colorBy = "cellColData",
  name = "MyeloidU"
)