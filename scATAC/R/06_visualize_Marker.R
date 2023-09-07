library(patchwork)
# load ArchR project
proj_path <- "Result/SIM12d" # 之前的 output dir
proj <- readRDS(paste0(proj_path, "/Save-ArchR-Project.rds"))
source("bin/umapForArchR.R")

p_umap <- umapForArchR(
  ArchR.proj = proj,
  embed.umap = "Tile_4",
  embed.cluster = "Tile_4",
  point.size = 1.5,
  title = "UMAP of SIM12d (Clustered by GeneScoreMatrix)"
)

p_umap_peaks <- umapForArchR(
  ArchR.proj = proj,
  embed.umap = "UMAP_Peaks",
  embed.cluster = "Cluster_Peaks",
  point.size = 1.5,
  title = "UMAP of SIM12d (Clustered by PeaksMatrix)"
)


pdf(paste0(proj_path, "/Plots/umap.pdf"), width = 14)
p_umap + p_umap_peaks
dev.off()


library(readxl)

marker_file <- "./Result/CIM22d/MarkerGene/Tile_4_MarkerGene.xlsx"

sheets <- excel_sheets(marker_file)
markerList <- purrr::map(
  sheets, ~ read_excel(marker_file, sheet = .)
)
names(markerList) <- sheets
tmp <- sapply(markerList, nrow)

markers <- do.call(rbind, markerList)

genes <- markers$name
gene_clu <- rep(names(tmp), tmp)

GSmat <- readRDS(paste0(proj_path, "/CIM22d_CY_GeneScoreMatrix.rds"))

source("bin/dotplotForArchR.R")

p_dot <- dotplotForArchR(
  geneScoreMatrix = GSmat, # 必要
  name = "Tile_4", # 必要
  # showclusters = c("C1", "C2", "C3", "C4", "C5"),
  genes = genes, # 必要
  rowsp = gene_clu,
  geneFontsize = 10,
  dotsize = 5,
  dotcolor = c("#2f0154", "#19b99f", "#f58114")
)

pdf(paste0(proj_path, "/Plots/Marker_dotplot_tile4.pdf"), height = 48)
p_dot
dev.off()
# proj <- addImputeWeights(
#   ArchRProj = proj,
#   reducedDims = "Tile_6"
#   # sampleCells = length(getCellNames(proj))
# )

markerGenes <- c("Os02g0809800")

# p <- plotEmbedding(
#   ArchRProj = proj,
#   colorBy = "GeneScoreMatrix",
#   name = markerGenes,
#   embedding = "UMAP_Tile_6",
#   rastr = FALSE,
#   imputeWeights = getImputeWeights(proj)
# )

p <- plotEmbedding(
  ArchRProj = proj,
  colorBy = "GeneScoreMatrix",
  name = markerGenes,
  embedding = "UMAP_Tile_6",
  rastr = FALSE,
  plotAs = "point",
  size = 1,
  pal = ArchRPalettes$whiteBlue,
  # baseSize = 20,
  quantCut = c(0.01, 0.99),
  imputeWeights = NULL
)


p2 <- lapply(p, function(x) {
  x + guides(color = "none", fill = "none") +
    theme_ArchR(baseSize = 18) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
})

pa <- do.call(cowplot::plot_grid, c(list(ncol = 2), p2))

pdf("c6pericycle.pdf", width = 15, height = 16)
pa
dev.off()

p3 <- plotBrowserTrack(
  ArchRProj = proj,
  groupBy = "Cluster_Tile_6",
  geneSymbol = markerGenes,
  upstream = 7000,
  downstream = 7000
)

grid::grid.newpage()
grid::grid.draw(p3$Os07g0684000)