Clu.sankey <- function(ArchR.proj = proj,
                       iter.Name1 = "Cluster_Tile_2",
                       iter.Name2 = "Cluster_Tile_4",
                       umap.Name1 = "UMAP_Tile_2",
                       umap.Name2 = "UMAP_Tile_4",
                       fontSize = 18,
                       height = NULL,
                       width = NULL) {
  if (!"networkD3" %in% installed.packages()[, 1]) {
    stop(
      "You need the packages: networkD3", "\n",
      "Try: install.packages(\"networkD3\")"
    )
  }
  if (!"scatterD3" %in% installed.packages()[, 1]) {
    stop(
      "You need the packages: scatterD3", "\n",
      "Try: install.packages(\"scatterD3\")"
    )
  }
  # if (!"manipulateWidget" %in% installed.packages()[, 1]) {
  #   stop(
  #     "You need the packages: manipulateWidget", "\n",
  #     "Try: install.packages(\"manipulateWidget\")"
  #   )
  # }

  cellColData <- ArchR.proj@cellColData
  embData <- ArchR.proj@embeddings

  if (!all(c(iter.Name1, iter.Name2) %in% names(cellColData))) {
    stop("? Please check the iter.Names: ", iter.Name1, "; ", iter.Name2)
  }
  if (!all(c(umap.Name1, umap.Name2) %in% names(embData))) {
    stop("? Please check the iter.Names: ", umap.Name1, "; ", umap.Name2)
  }

  col <- c(
    "#cc1616", "#d0595a", "#f786a8", "#701a12",
    "#faad89", "#f77c43", "#f85408", "#f01778",
    "#b7889d", "#ecd0d0", "#17e7ee", "#8dcaec",
    "#04a5fc", "#66c2a5", "#3bd608", "#a9e097",
    "#2c6917", "#cab9f5", "#df8633", "#05968e",
    "#5e6166", "#1872a3", "#7774c2", "#393b78",
    "#0b13f1", "#a00b98", "#63065d", "#2e012e"
  )

  its <- as.data.frame(
    cbind(
      cellColData[iter.Name1],
      cellColData[iter.Name2]
    )
  )

  cat(gsub("Length:", "Cell Numbers: ", summary(its)[1]), "\n")
  cat(iter.Name1, ":", unique(its[, 1]), "\n")
  cat(iter.Name2, ":", unique(its[, 2]), "\n")
  a <- as.data.frame(table(its[, 2]))[, 2]

  links <- data.frame(
    source = paste0(iter.Name1, "_", its[, 1]),
    target = paste0(iter.Name2, "_", its[, 2]),
    value = rep(a, a) / 100
  )

  nodes <- data.frame(
    name = c(
      as.character(links$source),
      as.character(links$target)
    ) %>% unique()
  )

  links$IDsource <- match(links$source, nodes$name) - 1
  links$IDtarget <- match(links$target, nodes$name) - 1


  p1 <- networkD3::sankeyNetwork(
    Links = links,
    Nodes = nodes,
    Source = "IDsource",
    Target = "IDtarget",
    Value = "value",
    NodeID = "name",
    fontSize = fontSize,
    height = height,
    width = width,
    sinksRight = FALSE
  )

  umap1 <- embData[[umap.Name1]]$df
  umap1 <- cbind(umap1, proj@cellColData[iter.Name1])
  colnames(umap1) <- c("UMAP_1", "UMAP_2", "Cluster")

  umap1_color <- col[sample(1:28)[seq.int(unique(umap1$Cluster))]]
  p2 <- scatterD3::scatterD3(
    data = umap1,
    x = UMAP_1, y = UMAP_2,
    col_var = Cluster,
    colors = umap1_color,
    point_opacity = 0.7
  )

  umap2 <- embData[[umap.Name2]]$df
  umap2 <- cbind(umap2, proj@cellColData[iter.Name2])
  colnames(umap2) <- c("UMAP_1", "UMAP_2", "Cluster")

  umap2_color <- col[sample(1:28)[seq.int(unique(umap2$Cluster))]]

  p3 <- scatterD3::scatterD3(
    data = umap2,
    x = UMAP_1, y = UMAP_2,
    col_var = Cluster,
    colors = umap2_color,
    point_opacity = 0.7
  )
  # p4 <- manipulateWidget::combineWidgets(p2, p3)
  return(list("sankey" = p1, "umap1" = p2, "umap2" = p3))
}
