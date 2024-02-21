if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# color
devtools::install_github("caleblareau/BuenColors")


# 安装cisTopic
devtools::install_github("aertslab/cisTopic")


# 安装ArchR
devtools::install_github("GreenleafLab/ArchR",
                         ref = "master",
                         repos = BiocManager::repositories())


library(ArchR)
ArchR::installExtraPackages()


install.packages("openxlsx")
install.packages("pheatmap")
install.packages("magick")

# For seurat 3.0
# packageurl <- "https://cran.r-project.org/src/contrib/Archive/Matrix/Matrix_1.3-2.tar.gz"
# install.packages(packageurl, repos=NULL, type="source")


# 准备archR依赖的注释


# 创建 genome annotation
if (!requireNamespace("BSgenome.Athaliana.TAIR.TAIR9", quietly = TRUE)) {
  BiocManager::install("BSgenome.Athaliana.TAIR.TAIR9", type="source")
  
}
library(BSgenome.Athaliana.TAIR.TAIR9)

blackList <- rtracklayer::import.bed("data/ref/blacklist.bed")

genomeAnnotation <- createGenomeAnnotation(genome = BSgenome.Athaliana.TAIR.TAIR9,
                                           filterChr = c("ChrM", "ChrC")
)
genomeAnnotation$blacklist <- blackList


# 创建gene annotation
# ./gff3_to_cellranger_gtf Araport11_GFF3_genes_transposons.201606.gff tmp.gtf
# uniq tmp.gtf > Araport11_GFF3_genes_transposons.201606.gtf
txdb <- GenomicFeatures::makeTxDbFromGFF('data/ref/Araport11_GFF3_genes_transposons.2016062.gtf')

# 一定要有symbol, 因为ArchR会在算gene score的时候过滤掉没有symbol的基因
genes <- genes(txdb)
genes$symbol <- genes$gene_id

exons <- unlist(exonsBy(txdb, "gene"))
mcols(exons)['gene_id'] <- names(exons)

transcripts <- transcripts(txdb)

tss <- GRanges(seqnames =seqnames(transcripts),
               ranges = IRanges(start(transcripts), width = 1),
               strand = strand(transcripts))
mcols(tss)['tx_id'] <- transcripts$tx_id
mcols(tss)['tx_name'] <- transcripts$tx_name


geneAnnotation <- createGeneAnnotation(genes = genes, exons=exons,TSS=tss)

# 保存
save(genomeAnnotation, geneAnnotation, file = "data/ref/TAIR10.Rdata")