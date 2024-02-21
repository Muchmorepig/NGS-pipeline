# install.packages("https://github.com/xuzhougeng/org.Osativa.eg.db/releases/download/v0.01/org.Osativa.eg.db.tar.gz",
#     repos = NULL,
#     type = "source"
# )
# library(org.Osativa.eg.db)
# orgdb <- org.Osativa.eg.db
library(clusterProfiler)
library(data.table)
library(magrittr)
library(readxl)
library(fs)
library(org.At.tair.db)
orgdb <- org.At.tair.db

os_at <- fread(
    "/data3/wanglab/wmc/reference/Oryza_sativa/annotation/IRGSP_OsAt_Orthologous.tsv"
)
marker_file <- "./Result/SIM12d/Marker/MarkerGene_Tile4.xlsx"

sheets <- excel_sheets(marker_file)
markerList <- purrr::map(
    sheets, ~ read_excel(marker_file, sheet = .)
)
names(markerList) <- sheets
sapply(markerList, nrow)

goall <- lapply(sheets, function(x) {
    if (nrow(markerList[[x]]) == 0) {
        cat("Cluster:", x, " 没有Marker\n")
        return(NULL)
    }
    gene <- unique(
        os_at[RAP_ID %in% markerList[[x]]$name, TAIR]
    )
    if (length(gene) == 0) {
        cat("Cluster:", x, " 没有拟南芥同源\n")
        return(NULL)
    }
    cat("enrichGO for Cluster:", x, "\n")
    gores <- enrichGO(
        gene = gene,
        OrgDb = orgdb,
        keyType = "TAIR",
        ont = "BP",
        readable = TRUE,
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2
    ) %>%
        clusterProfiler::simplify(
            .,
            cutoff = 0.7,
            by = "p.adjust",
            select_fun = min
        ) %>%
        as.data.frame()
    if (nrow(gores) == 0) {
        NULL
    } else {
        gores
    }
})
names(goall) <- sheets
# 保存
openxlsx::write.xlsx(goall, paste0(path_dir(marker_file), "/Tile4_MarkerGene_go.xlsx"))


sel <- names(unlist(sapply(goall, nrow)))

# 可视化
source("bin/goplot.R")

gobarplot(goall[[1]])
gobarplot(goall[[2]])
gobarplot(goall[[3]])

gonetplort(goall[[1]],
    showCategory = 5,
    circular = TRUE,
    colorEdge = TRUE,
    node_label = "gene",
    alias = NULL
)

p <- lapply(sel, function(x) {
    gobarplot(goall[[x]], tt = x, size = 4)
})

pall <- do.call(cowplot::plot_grid, c(list(ncol = 2), p))


pdf("Result/SIM12d/Plots/markerGene_go.pdf", height = 28, width = 12)
pall
dev.off()
