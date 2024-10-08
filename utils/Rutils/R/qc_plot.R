#' @import ggplot2
#' @import dplyr
#' @import patchwork
#' @export

qc_plot <- function(obj = scrna,
                    max_nFeature = 6000,
                    min_nFeature = 500,
                    max_nCount = 100000,
                    min_nCount = 1000,
                    color = c("percent.mito", "percent.chlo")) {
  # Helper function to create violin plot with thresholds
  create_violin_plot <- function(obj, feature, min_val, max_val, add_thresholds = TRUE) {
    p <- suppressMessages(VlnPlot(obj, features = feature, layer = "counts")) +
      scale_y_continuous(expand = c(0, 0)) +
      theme(
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none"
      )
    if (add_thresholds) {
      p <- p +
        geom_hline(yintercept = min_val, linetype = "dashed") +
        geom_hline(yintercept = max_val, linetype = "dashed")
    }

    return(p)
  }

  p_feature <- create_violin_plot(obj, "nFeature_RNA", min_nFeature, max_nFeature)
  p_count <- create_violin_plot(obj, "nCount_RNA", min_nCount, max_nCount)

  p_scatter <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
    theme(legend.position = "none") +
    geom_hline(yintercept = c(min_nFeature, max_nFeature), linetype = "dashed") +
    geom_vline(xintercept = c(min_nCount, max_nCount), linetype = "dashed")

  p_main <- p_feature | p_count | p_scatter

  if (!is.null(color) && all(color != "") && length(color) > 0) {
    plist_color <- lapply(color, function(f) create_violin_plot(obj, f, min_nFeature, max_nFeature, add_thresholds = FALSE))

    dd <- as_tibble(obj[[]], rownames = "Cell.Barcode")

    plist2 <- lapply(color, function(col) {
      dd <- dd %>% arrange(!!sym(col))
      ggplot(dd, aes(x = nCount_RNA, y = nFeature_RNA, color = !!sym(col))) +
        geom_point(size = 1, alpha = 0.9) +
        theme_classic() +
        scale_color_gradientn(
          colors = c("grey", "blue", "green2", "red", "yellow"),
          guide = guide_legend(title = col, title.position = "top")
        ) +
        scale_y_continuous(expand = c(0, min(min_nFeature - 200, 0.05 * max_nFeature))) +
        theme(legend.position = c(0.8, 0.3))
    })

    p_color <- wrap_plots(plist_color) + plot_layout(ncol = length(color))
    p_color <- p_color & NoLegend() & theme(axis.title = element_blank(), axis.text.x = element_blank())
    p_gradient <- wrap_plots(plist2) + plot_layout(ncol = length(color))

    return(p_main / (p_color | p_gradient))
  }

  return(p_main)
}
