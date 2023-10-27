require(ggplot2)
require(patchwork)
require(magrittr)

qc_plot <- function(obj = scrna,
                    max_nFeature = 6000,
                    min_nFeature = 500,
                    feature1 = "nCount_RNA",
                    feature2 = "nFeature_RNA",
                    color = c("percent.mito", "percent.chlo")) {
  p1 <- suppressMessages(
    Seurat::VlnPlot(obj,
      features = c("nFeature_RNA")
    ) &
      geom_hline(yintercept = min_nFeature, linetype = "dashed") &
      geom_hline(yintercept = max_nFeature, linetype = "dashed")
  )

  plist1 <- lapply(c(feature1, color), function(f) {
    suppressMessages(
      Seurat::VlnPlot(obj, features = f)
    )
  })

  p2 <- suppressMessages(
    (p1 | plist1[[1]] | plist1[[2]] | plist1[[3]]) &
      scale_y_continuous(expand = c(0, 0)) &
      theme(
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none"
      )
  )


  p3 <- Seurat::FeatureScatter(obj,
    feature1 = feature1,
    feature2 = feature2
  ) +
    ggplot2::theme(legend.position = "none")

  dd <- dplyr::as_tibble(
    obj[[]],
    rownames = "Cell.Barcode"
  )

  plist2 <- lapply(color, function(col) {
    dd <- dd[as.data.frame(dd) %>%
      .[, col] %>%
      order(), ]
    dd %>%
      ggplot(aes(nCount_RNA, nFeature_RNA, colour = unlist(as.vector(dd[, col])))) +
      geom_point(size = 1, alpha = 0.9) +
      theme_classic() +
      scale_color_gradientn(
        colors = c("grey", "blue", "green2", "red", "yellow"),
        guide = guide_legend(
          title = col,
          title.position = "top"
          # direction = "horizontal"
        )
      ) +
      geom_hline(yintercept = min_nFeature, linetype = "dashed") +
      geom_hline(yintercept = max_nFeature, linetype = "dashed") +
      scale_y_continuous(expand = c(0, min(min_nFeature - 200, 0.05 * max_nFeature))) +
      theme(legend.position = c(0.8, 0.3))
  })

  p4 <- p3 + plist2[[1]] + plist2[[2]]
  p2 / p4
}
