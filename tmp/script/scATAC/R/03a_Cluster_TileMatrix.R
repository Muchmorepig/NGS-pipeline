suppressMessages(library(ArchR))
library(parallel)
library(magrittr)
library(fs)
addArchRThreads(threads = 32)

set.seed(1234) # 可加可不加
# load ArchR project
proj_path <- "Result/SIM16d" # 之前的 output dir
proj <- loadArchRProject(path = proj_path)

# 设置不同的迭代次数并命名，然后会按照不同的次数迭代降维, 可能有点耗时 ------
iter_vector <- c(2, 4, 6)
names(iter_vector) <- c("Tile_2", "Tile_4", "Tile_6")

for (x in seq_along(iter_vector)) {
    iter <- iter_vector[x]
    lsi_name <- names(iter_vector[x])
    cat("iterations: ", iter, "\t", "name: ", lsi_name, "\n")
    proj %<>%
        addIterativeLSI(
            ArchRProj = .,
            useMatrix = "TileMatrix",
            name = lsi_name,
            iterations = iter,
            clusterParams = list(
                resolution = c(0.8),
                sampleCells = 10000,
                n.start = 10
            ),
            force = TRUE
        ) %>%
        addClusters(
            input = .,
            reducedDims = lsi_name,
            # name = paste0("Cluster_", lsi_name),
            name = lsi_name,
            force = TRUE
        ) %>%
        addUMAP(
            ArchRProj = .,
            reducedDims = lsi_name,
            nNeighbors = 40,
            minDist = 0.4,
            metric = "cosine",
            # name = paste0("UMAP_", lsi_name),
            name = lsi_name,
            force = TRUE
        )
    # %>%
    # addTSNE(
    #     ArchRProj = .,
    #     reducedDims = lsi_name,
    #     name = paste0("TSNE_", lsi_name),
    #     perplexity = 30
    # )
}

# save ArchRProject -------
saveArchRProject(proj)


source("./bin/umapForArchR.R")

umapForArchR(ArchR.proj = proj, embed.umap = "res_08", embed.cluster = "res_08")

## 可以按照Cluster做一些基本统计
library(dplyr)
as_tibble(proj@cellColData) %>%
    group_by(Tile_6) %>%
    summarise(min(TSSEnrichment), mean(TSSEnrichment), median(TSSEnrichment), max(TSSEnrichment)) %>%
    ungroup()


# 画个 Sankey Plot，看看不同的迭代次数（没改resolution）的影响
source("~/script/tips_func/sankey_umap_ForArchR.R")

p1 <- Clu.sankey_ummap(proj,
    iter.Name1 = "Tile_2",
    iter.Name2 = "Tile_4",
    umap.Name1 = "Tile_2",
    umap.Name2 = "Tile_4"
)

p2 <- Clu.sankey_ummap(proj,
    umap.Name1 = "Tile_4",
    umap.Name2 = "Tile_6",
    iter.Name1 = "Tile_4",
    iter.Name2 = "Tile_6"
)

# library(manipulateWidget)
# p3 <- combineWidgets(p1$sankey, p2$sankey,
#     ncol = 2,
#     title = "<h1 align='center'>SankeyPlot: Comparing Cluster used different iterations</h1>",
#     header = "<h3 align='center' style='background-color:#cdcdcd;'>
# 左: 迭代2次与4次; 右: 迭代4次与8次."
# )

# dir_create(paste0(proj_path, "/Plots"))
# 保存 html
htmlwidgets::saveWidget(p1$sankey,
    file = "Plots/sankey_iter24.html",
    selfcontained = FALSE
)
htmlwidgets::saveWidget(p2$sankey,
    file = "Plots/sankey_iter48.html",
    selfcontained = FALSE
)

# htmlwidgets::saveWidget(p4,
#     file = paste0(proj_path, "/Plots/umap.html"),
#     selfcontained = FALSE
# )

htmlwidgets::saveWidget(p1$umap1,
    file = "Plots/umap_iter2.html",
    selfcontained = FALSE
)

htmlwidgets::saveWidget(p1$umap2,
    file = "Plots/umap_iter4.html",
    selfcontained = FALSE
)

htmlwidgets::saveWidget(p2$umap2,
    file = "Plots/umap_iter8.html",
    selfcontained = FALSE
)

# you save it as an html
# saveNetwork(sn, "sn.html")
# library(webshot)
# # you convert it as png
# webshot("sn.html", "sn.png", vwidth = 1000, vheight = 900)


# plotly visualization
# library(plotly)
# df <- proj@embeddings$TSNE$df
# colnames(df) <- c("TSNE_1", "TSNE_2")
# df <- proj@embeddings$UMAP$df
# colnames(df) <- c("UMAP_1", "UMAP_2")
# df$group <- proj$Clusters

# library(scatterD3)
# emb <- proj@embeddings
# iter1 <- emb[["UMAP_Tile_4"]]$df
# iter1 <- cbind(iter1, proj@cellColData["Cluster_Tile_4"])
# colnames(iter1) <- c("UMAP_1", "UMAP_2", "Cluster")


# col <- c(
#     "#d72828", "#d0595a", "#f07c9f", "#c2534b", "#faad89", "#f57338", "#f85408",
#     "#f01778", "#b7889d", "#ecd0d0", "#17e7ee", "#8dcaec", "#04a5fc", "#66c2a5",
#     "#3bd608", "#a9e097", "#2c6917", "#d5c5ff", "#df8633", "#05968e", "#5e6166",
#     "#1872a3", "#7774c2", "#393b78", "#0b13f1", "#a00b98", "#460a42", "#070007"
# )

# scatterD3(
#     data = iter1,
#     x = UMAP_1, y = UMAP_2,
#     col_var = Cluster,
#     colors = col[seq.int(unique(iter1$Cluster))],
#     point_opacity = 0.7
# )