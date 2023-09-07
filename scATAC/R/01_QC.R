suppressMessages(library(ArchR))
library(parallel)
library(fs)
addArchRThreads(threads = 32)
addArchRChrPrefix(chrPrefix = FALSE)
load("MSU7.Rdata")

# 数据位置，样本名称，输出目录 --------
inputFiles <- "rawdata/Os5dR_MX_fragments.tsv.gz"
sampleNames <- "Os5dR"

output_dir <- "temp"
if (!dir_exists(output_dir)) dir_create(output_dir)

# pre-filetration --------
arrowFiles <- createArrowFiles(
    inputFiles = inputFiles,
    sampleNames = sampleNames,
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

proj <- ArchRProject(
    ArrowFiles = arrowFiles,
    outputDirectory = output_dir,
    geneAnnotation = geneAnnotation,
    genomeAnnotation = genomeAnnotation,
    copyArrows = FALSE
)

df <- getCellColData(proj, select = c("log10(nFrags)", "TSSEnrichment", "Sample"))
p <- ggPoint(
    x = df[, 1],
    y = df[, 2],
    colorDensity = TRUE,
    continuousSet = "sambaNight",
    size = 0.75,
    xlabel = "Log10 Unique Fragments",
    ylabel = "TSS Enrichment",
    xlim = c(log10(100), log10(100000)),
    ylim = c(0, 4),
    rastr = TRUE
) + geom_hline(yintercept = 1.5, lty = "dashed") + geom_vline(xintercept = log10(1000), lty = "dashed")


plotPDF(p, name = paste0(sampleNames, "_qc"))

# 删除不必要的文件
file_delete(arrowFiles)
dir_delete(output_dir)
