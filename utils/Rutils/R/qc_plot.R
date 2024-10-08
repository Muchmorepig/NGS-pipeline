#' @export

require(ggplot2)
require(patchwork)
require(dplyr)
qc_plot <- function(obj = scrna,
                    max_nFeature = 6000,
                    min_nFeature = 500,
                    max_nCount = 100000,
                    min_nCount = 1000,
                    color = c("percent.mito", "percent.chlo")) {
  # Helper function to create violin plot with thresholds
  create_violin_plot <- function(obj, feature, min_val, max_val, add_thresholds = TRUE) {
    p <- suppressMessages(VlnPlot(obj, features = feature, layer = "counts"))

    if (add_thresholds) {
      p <- p +
        geom_hline(yintercept = min_val, linetype = "dashed") +
        geom_hline(yintercept = max_val, linetype = "dashed")
    }

    return(p)
  }

  p_feature <- create_violin_plot(obj, "nFeature_RNA", min_nFeature, max_nFeature)
  p_count <- create_violin_plot(obj, "nCount_RNA", min_nCount, max_nCount)

  p3 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
    theme(legend.position = "none")

  p2 <- (p_feature | p_count) &
    scale_y_continuous(expand = c(0, 0)) &
    theme(
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.line.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "none"
    )
  p1 <- p2 | (p3 + geom_hline(yintercept = min_nFeature, linetype = "dashed") +
    geom_hline(yintercept = max_nFeature, linetype = "dashed") + geom_vline(xintercept = min_nCount, linetype = "dashed") +
    geom_vline(xintercept = max_nCount, linetype = "dashed"))

  if (!is.null(color) && color != "") {
    plist_color <- lapply(color, function(f) create_violin_plot(obj, f, min_nFeature, max_nFeature, add_thresholds = FALSE))

    dd <- as_tibble(obj[[]], rownames = "Cell.Barcode")

    plist2 <- lapply(color, function(col) {
      dd <- dd %>% arrange(!!sym(col))
      ggplot(dd, aes_string(x = "nCount_RNA", y = "nFeature_RNA", color = col)) +
        geom_point(size = 1, alpha = 0.9) +
        theme_classic() +
        scale_color_gradientn(
          colors = c("grey", "blue", "green2", "red", "yellow"),
          guide = guide_legend(title = col, title.position = "top")
        ) +
        geom_hline(yintercept = min_nFeature, linetype = "dashed") +
        geom_hline(yintercept = max_nFeature, linetype = "dashed") +
        scale_y_continuous(expand = c(0, min(min_nFeature - 200, 0.05 * max_nFeature))) +
        theme(legend.position = c(0.8, 0.3))
    })

    p_color <- wrap_plots(plist_color) + plot_layout(ncol = length(color))
    p_color <- p_color & NoLegend() & theme(axis.title = element_blank(), axis.text.x = element_blank())
    p_gradient <- wrap_plots(plist2) + plot_layout(ncol = length(color))

    return(p1 / (p_color | p_gradient))
  }

  return(p1)
}

