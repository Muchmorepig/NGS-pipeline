# library(parallel)
# library(ArchR)
# library(Seurat)
# genes_bc <- read.table("./rawdata/CIM22d_genes_bc.bed",
#   sep = "\t",
#   as.is = c(4, 7),
#   header = FALSE
# )
#
# bc <- genes_bc[, 7]
# bc_split <- strsplit(as.character(bc), ";")
# bc_split_vec <- unlist(bc_split)
# bc_unique <- unique(bc_split_vec)
# bc_counts <- table(bc_split_vec)
# bc_filt <- names(bc_counts)[bc_counts > 800]
# barcodes <- bc_filt
#
# gene_counts <- makeFeatureMatrix(genes_bc, barcodes)
# gene_counts <- gene_counts[order(rownames(gene_counts)), ]
# colnames(gene_counts) <- paste0("atac_", colnames(gene_counts))
#
# scrna <- read10X(
#   sample.dirs = list("/data5/cy_data/scRNA_Seq/202111_SC20210915S1_OSCIM22D_EXTRA/run_cellranger_count/run_count_OSCIM22DEXTRA"),
#   sample.names = list("rna")
# )
# dat <- list(atac = gene_counts, rna = scrna)

# proj <- "./Result/Merge"
# proj <- ArchR::loadArchRProject(path = proj)
# gs_df <- ArchR::getMatrixFromProject(proj)
#
# saveRDS(gs_df, "./Result/Merge/merge_CY_GeneScoreMatrix.rds")

#
# pdf("Gene_Loadings.pdf")
# plotGeneLoadings(int_dat, return.plots = FALSE)
# dev.off()


# p <- plotGene(int_dat,
#   "Os10g0542100",
#   axis.labels = c("UMAP 1", "UMAP 2"),
#   cols.use = c("white", "red"),
#   pt.size = 0.3,
#   return.plots = TRUE
# )
#
# p$atac | p$rna
#
# b <- a[[1]]$patches$plots
# b <- b[[2]]
# grep("Os10g0542", rownames(b$patches$plots[[2]]$data), value = T)


library(magrittr)
library(rliger)
gs_df <- readRDS("./Result/CIM22d/CIM22d_CY_GeneScoreMatrix.rds")
# gs_df <- readRDS("./Result/Os5dR/Os5dR_GeneScoreMatrix.rds")
# gs_df <- readRDS("./Result/Merge/merge_CY_GeneScoreMatrix.rds")
scatac <- assay(gs_df)
row.names(scatac) <- rowData(gs_df)$name

# scrna <- Seurat::Read10X(sc_rna_path)
# scrna_path <- "./rawdata/osRoots.Rds"
scrna_path <- "/data5/cy_data/scRNA_Seq/202111_SC20210915S1_OSCIM22D_EXTRA/OSCIM22D_EXTRA_downstream/cim22d_extra_f8k_hvg1.rds"
# scrna_path <- "/data5/cy_data/scRNA_Seq/202203_oscim22dsim4dsimm12d_integrate/cim22d_sim4d_sim12d.rds"
scrna <- readRDS(scrna_path)
scrna <- scrna@assays$RNA@counts

dat <- list(atac = scatac, rna = scrna)

int_dat <- createLiger(dat)
# take.gene.union = T

rm(gs_df, scatac, scrna_path, scrna, dat)
gc()


int_dat %<>% rliger::normalize() %>%
  selectGenes(
    .,
    datasets.use = 2,
    var.thresh = 0.1
  )

length(int_dat@var.genes)

int_dat %<>% scaleNotCenter() %>% optimizeALS(., k = 20)

int_dat %<>% quantile_norm() %>%
  louvainCluster(., resolution = 0.5) %>%
  runUMAP(
    .,
    distance = "cosine",
    n_neighbors = 30,
    min_dist = 0.3
  )

p <- plotByDatasetAndCluster(
  int_dat,
  pt.size = 0.6,
  axis.labels = c("UMAP 1", "UMAP 2"),
  return.plots = TRUE
)

pdf("CIM22d.pdf", height = 8, width = 16)
p[[1]] | p[[2]]
dev.off()

pp <- function(gene_chr, liger_dat = int_dat, pdf = "tmp.pdf") {
  p_list <- lapply(gene_chr, function(x) {
    p <- plotGene(liger_dat,
      x,
      axis.labels = c("UMAP 1", "UMAP 2"),
      cols.use = c("white", "purple"),
      pt.size = 0.3,
      return.plots = TRUE
    )
    return(p)
  })

  names(p_list) <- gene_chr

  tmp <- lapply(gene_chr, function(x) {
    plot_grid(p_list[[x]]$atac, p_list[[x]]$rna, ncol = 2)
  })

  pdf(pdf, width = 16)
  print(tmp)
  dev.off()
}

gene_chr <- c(
  "Os10g0122600", "Os03g0150800", "Os10g0578200", "Os06g0108600",
  "Os06g0184000", "Os01g0248900", "Os02g0595900", "Os12g0637100",
  "Os04g0452700", "Os03g0107300", "Os03g0279200", "Os05g0438700",
  "Os08g0490900", "Os02g0805200", "Os02g0829100", "Os08g0512600",
  "Os01g0805600", "Os12g0555600", "Os04g0486500", "Os02g0810200",
  "Os06g0127800", "Os07g0105700", "Os03g0820500", "Os04g0540900",
  "Os04g0546800", "Os10g0467800", "Os09g0422500", "Os04g0536500",
  "Os03g0640800", "Os06g0614000", "Os06g0676000", "Os01g0236300",
  "Os01g0878700", "Os01g0847300", "Os03g0216700", "Os04g0445000"
)

"OsCESA9 Os09g0422500 vc"
"OsSKOR Os04g0445000 pe"
"OsNAR2.1 Os02g0595900 ep"

pp(gene_chr = gene_chr, liger_dat = int_dat, pdf = "tmp.pdf")


# 用liger的找一些marker
library(dplyr)

int_dat_wilcoxon <- runWilcoxon(
  int_dat,
  data.use = "all",
  compare.method = "clusters"
)

int_dat_wilcoxon %<>% subset(padj < 0.05 & logFC > 4)

cluster_df <- data.frame(cluster = int_dat@clusters)
cluster_df$dataset <- int_dat@cell.data$dataset

cc <- unique(cluster_df$cluster) %>% levels()

genes <- lapply(cc, function(x) {
  markers <- int_dat_wilcoxon[int_dat_wilcoxon$group == x, ]
  markers <- markers[order(markers$padj), ]
  markers <- markers[1:2, ]
  return(markers$feature)
}) %>% unlist()


pp(gene_chr = genes, liger_dat = int_dat, pdf = "tmp_cim22d.pdf")



# 获取atac-seq和rna-seq的聚类标记 在整合的结果里面展示

## 获取rna-seq的聚类信息
seu.obj <- readRDS(scrna_path)

rna_cluster <- seu.obj$seurat_clusters
rna_cluster_df <- as.data.frame(rna_cluster)
rna_cluster_df$barcode <- row.names(rna_cluster_df)

tsne_df <- as.data.frame(int_dat@tsne.coords)
tsne_df$barcode <- row.names(tsne_df)
tsne_df <- merge(tsne_df, rna_cluster_df, by = "barcode")

library(ggplot2)
ppp <- ggplot(tsne_df, aes(x = V1, y = V2)) +
  geom_point(aes(color = rna_cluster), size = 0.7) +
  scale_color_manual(values = BuenColors::jdb_palette("corona")) +
  theme_classic() +
  labs(x = "UMAP 1", y = "UMAP 2")

ppp

## 获取atac_cluster的cluster信息
proj <- ArchR::loadArchRProject("./Result/CIM22d")
atac_cluster <- proj@cellColData$Tile_6
names(atac_cluster) <- row.names(proj@cellColData)
atac_cluster_df <- as.data.frame(atac_cluster)
atac_cluster_df$barcode <- row.names(atac_cluster_df)

tsne_df <- as.data.frame(int_dat@tsne.coords)
tsne_df$barcode <- row.names(tsne_df)

tsne_df <- merge(tsne_df, atac_cluster_df, by = "barcode")

ppp1 <- ggplot(tsne_df, aes(x = V1, y = V2)) +
  geom_point(aes(color = atac_cluster), size = 0.7) +
  scale_color_manual(values = BuenColors::jdb_palette("corona")) +
  theme_classic() +
  labs(x = "UMAP 1", y = "UMAP 2")



pdf("CIM22d_pri.pdf", height = 8, width = 16)
ppp | ppp1
dev.off()

#
# df_cellranger <- data.table::fread("../atac2/outs/singlecell.csv")
# cellranger_cell <- df_cellranger$barcode[df_cellranger$is__cell_barcode == 1]
#
#
# mycell <- substring(proj$cellNames, 8)
#
# gplots::venn(list(ArchR = cellranger_cell, CellRanger = mycell))
# g <- venn.diagram(
#   x = list(cellranger_cell, mycell),
#   category.names = c("ArchR", "CellRanger"),
#   filename = NULL,
#   output = TRUE
# )
# grid.newpage()
# grid.draw(g)

df_glue <- read.csv("./scglue/scglue_combined.csv")
colnames(df_glue)[1] <- "barcode"
df_glue$leiden <- as.character(df_glue$leiden)
tsne_df <- as.data.frame(int_dat@tsne.coords)
tsne_df$barcode <- gsub("CIM22d#", "", row.names(tsne_df))

table(tsne_df$barcode %in% df_glue$barcode)


ts <- merge(tsne_df, df_glue, by = "barcode")


ppp2 <- ggplot(ts, aes(x = V1, y = V2)) +
  geom_point(aes(color = leiden), size = 0.7) +
  scale_color_manual(values = BuenColors::jdb_palette("corona")) +
  theme_classic() +
  labs(x = "UMAP 1", y = "UMAP 2")

pdf("tt.pdf")
ppp2
dev.off()



glue <- read.csv("./scglue/scglue_combined_umap.csv")
glue$leiden <- as.character(glue$leiden)

ggplot(glue, aes(x = umap1, y = umap2)) +
  geom_point(aes(color = leiden), size = 0.7) +
  scale_color_manual(values = BuenColors::jdb_palette("corona")) +
  theme_classic() +
  labs(x = "UMAP 1", y = "UMAP 2")