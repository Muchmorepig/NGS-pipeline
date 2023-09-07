gs_df <- ArchR::getMatrixFromProject(proj)

gs_matrix <- assay(gs_df)
row.names(gs_matrix) <- rowData(gs_df)$name


gex_matrix <- Read10X(data.dir = "sim6d_rna/filtered_feature_bc_matrix/")


sim6d.data <- list(atac = gs_matrix, rna = gex_matrix)

int.sim6d <- createLiger(sim6d.data)
rm(gs_df, gs_matrix, gex_matrix, sim6d.data)
gc()


int.sim6d <- normalize(int.sim6d)
int.sim6d <- selectGenes(int.sim6d, datasets.use = 2)
int.sim6d <- scaleNotCenter(int.sim6d)

int.sim6d <- optimizeALS(int.sim6d, k = 20)

int.sim6d <- quantile_norm(int.sim6d)
int.sim6d <- louvainCluster(int.sim6d, resolution = 0.2)


int.sim6d <- runUMAP(int.sim6d, distance = "cosine", n_neighbors = 30, min_dist = 0.3)


p <- plotByDatasetAndCluster(int.sim6d, axis.labels = c("UMAP 1", "UMAP 2"), return.plots = TRUE)
p[[1]] | p[[2]]

p <- plotGene(int.sim6d, "AT2G17950", axis.labels = c("UMAP 1", "UMAP 2"), return.plots = TRUE)
p <- plotGene(int.sim6d, "AT1G64700", axis.labels = c("UMAP 1", "UMAP 2"), return.plots = TRUE)

p$atac | p$rna


cluster_df <- data.frame(cluster = int.sim6d@clusters)
cluster_df$dataset <- int.sim6d@cell.data$dataset



# 获取atac-seq和rna-seq的聚类标记 在整合的结果里面展示

## 获取rna-seq的聚类信息

rna_cluster <- seu.obj$seurat_clusters
rna_cluster_df <- as.data.frame(rna_cluster)
rna_cluster_df$barcode <- row.names(rna_cluster_df)

tsne_df <- as.data.frame(int.sim6d@tsne.coords)
tsne_df$barcode <- row.names(tsne_df)

tsne_df <- merge(tsne_df, rna_cluster_df, by = "barcode")


pp <- ggplot(tsne_df, aes(x = V1, y = V2)) +
  geom_point(aes(color = rna_cluster), size = 0.7) +
  scale_color_manual(values = BuenColors::jdb_palette("corona")) +
  theme_classic()

pp

pp1 <- pp + labs(x = "UMAP 1", y = "UMAP 2")
pp1

## 获取atac_cluster的cluster信息
atac_cluster <- proj@cellColData$Clusters
names(atac_cluster) <- row.names(proj@cellColData)
atac_cluster_df <- as.data.frame(atac_cluster)
atac_cluster_df$barcode <- row.names(atac_cluster_df)
## 获取
tsne_df <- as.data.frame(int.sim6d@tsne.coords)
tsne_df$barcode <- row.names(tsne_df)


tsne_df <- merge(tsne_df, atac_cluster_df, by = "barcode")

pp <- ggplot(tsne_df, aes(x = V1, y = V2)) +
  geom_point(aes(color = atac_cluster), size = 0.7) +
  scale_color_manual(values = BuenColors::jdb_palette("corona")) +
  theme_classic()

pp

# pp <- pp + geom_label(aes(x=mean_x,y=mean_y, label=color), text_anno)
pp2 <- pp + labs(x = "UMAP 1", y = "UMAP 2")

pp2

pp1 | pp2

#
df_cellranger <- data.table::fread("../atac2/outs/singlecell.csv")
cellranger_cell <- df_cellranger$barcode[df_cellranger$is__cell_barcode == 1]


mycell <- substring(proj$cellNames, 8)

gplots::venn(list(ArchR = cellranger_cell, CellRanger = mycell))
g <- venn.diagram(
  x = list(cellranger_cell, mycell),
  category.names = c("ArchR", "CellRanger"),
  filename = NULL,
  output = TRUE
)
grid.newpage()
grid.draw(g)