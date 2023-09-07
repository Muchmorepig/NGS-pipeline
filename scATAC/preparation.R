library(ArchR)
# ArchR::installExtraPackages()
library(TxDb.Osativa.MSU.xzg)
library(BSgenome.Osativa.MSU.xzg)

# 创建 genome annotation
# blackList <- rtracklayer::import.bed("ref/blacklist.bed")
# blackList <- NULL
genomeAnnotation <- createGenomeAnnotation(
    genome = BSgenome.Osativa.MSU.xzg,
    filterChr = c("ChrM", "ChrC")
)
# genomeAnnotation$blacklist <- blackList


# 创建gene annotation
txdb <- TxDb.Osativa.MSU.xzg

# 一定要有symbol, 因为ArchR会在算gene score的时候过滤掉没有symbol的基因
genes <- genes(txdb)
genes$symbol <- genes$gene_id

exons <- unlist(exonsBy(txdb, "gene"))
mcols(exons)["gene_id"] <- names(exons)

transcripts <- transcripts(txdb)

tss <- GRanges(
    seqnames = seqnames(transcripts),
    ranges = IRanges(start(transcripts), width = 1),
    strand = strand(transcripts)
)
mcols(tss)["tx_id"] <- transcripts$tx_id
mcols(tss)["tx_name"] <- transcripts$tx_name


geneAnnotation <- createGeneAnnotation(genes = genes, exons = exons, TSS = tss)

# 保存
save(genomeAnnotation, geneAnnotation, file = "MSU7.Rdata")