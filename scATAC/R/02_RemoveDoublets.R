suppressMessages(library(ArchR))
library(parallel)
library(magrittr)
library(fs)
addArchRThreads(threads = 32)
addArchRChrPrefix(chrPrefix = FALSE)

load("MSU7.Rdata")

# 数据位置，样本名称，输出目录 --------
inputFiles <- "rawdata/Os5dR_MX_fragments.tsv.gz"
if (!file_exists(inputFiles)) {
    stop("检查文件是否存在")
}

sampleNames <- "Os5dR"
output_dir <- paste0("Result/", sampleNames)
if (!dir_exists(output_dir)) dir_create(output_dir)

# 根据01_QC的图，调整参数，创建arrow文件，移动到 arrowFiles 目录 -------
minTSS <- 1.5
# minTSS <- 1.7
minFrags <- 10000 #3.75 10^3.75
# minFrags <- 6310 # 3.8
maxFrags <- 1e6

arrowFiles <- createArrowFiles(
    inputFiles = inputFiles,
    sampleNames = sampleNames,
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
    force = TRUE,
    cleanTmp = TRUE
)


doubScores <- addDoubletScores(
    input = arrowFiles,
    k = 10,
    knnMethod = "UMAP",
    LSIMethod = 1
)


proj <- ArchRProject(
    ArrowFiles = arrowFiles,
    outputDirectory = output_dir,
    copyArrows = TRUE, # This is recommened so that you maintain an unaltered copy for later usage.,
    geneAnnotation = geneAnnotation,
    genomeAnnotation = genomeAnnotation,
)

# Remove doublets
proj <- filterDoublets(proj)

# 保存
saveArchRProject(proj)
# 保存 GeneScoreMatrix
GSmat <- getMatrixFromProject(proj, "GeneScoreMatrix")
saveRDS(GSmat, paste0(output_dir, "/Os5dR_GeneScoreMatrix.rds"))

# ArchRProject 复制了一份 arrow 文件，把之前的删掉，节省地方 -------
file_size(arrowFiles)
file_delete(arrowFiles)
