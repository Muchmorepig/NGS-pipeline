if (!"dplyr" %in% (.packages())) library(dplyr, quietly = TRUE)
if (!"ggplot2" %in% (.packages())) library(ggplot2, quietly = TRUE)
if (!"gtools" %in% installed.packages()[, 1]) install.packages("gtools")

#' umap plot for Archr obj
#'
#' @param ArchR.proj ArchR project
#' @param mtx SummarizedExperiment project of PeakMatrix or GeneScoreMatrix
#' @param colorBy it can be one gene or Clusters Name in cellColData
#' @param embed.umap the shape of your plot, it can be obtained in getEmbedding
#' @param multiple bool, whether return a single plot with all Clusters or faced plot with single Cluster

umapForArchR <- function(ArchR.proj = proj,
                         embed.umap = "UMAP_Tile_8",
                         colorBy = "Cluster_Tile_8",
                         mtx = NULL,
                         marker.color = c(low = "#eafff2", mid = "#fad8bf", high = "#ff0000"),
                         palette = paletteer::paletteer_d("ggsci::category20c_d3"),
                         multiple = FALSE,
                         label.size = 3.6,
                         label.shadow = 0.25,
                         label.color = "#000103",
                         shadow.color = "#f5f0f0",
                         point.size = 0.8,
                         point.alpha = 0.8,
                         theme = c("mini", "clear"),
                         title = NULL,
                         title.size = 13,
                         legend.size = 10,
                         log2Norm = TRUE) {
    if (!embed.umap %in% names(ArchR.proj@embeddings)) stop("`embed.umap` must be inclued in `names(ArchR.proj@embeddings`")

    df <- ArchR.proj@embeddings[[embed.umap]]$df
    colnames(df) <- c("UMAP_Dim1", "UMAP_Dim2")
    # && is.null(marker)
    if (is.null(mtx)) {
        celcol <- colnames(ArchR.proj@cellColData)
        colorBy <- celcol[match(colorBy, celcol)]
        if (is.na(colorBy)) {
            stop("`colorBy` must be inclued in `colnames(ArchR.proj@cellColData)`")
        }
        df$Clusters <- as.vector(ArchR.proj@cellColData[, colorBy])
        df$Clusters <- factor(df$Clusters,
            levels = gtools::mixedsort(unique(df$Clusters))
        )
        if (multiple) {
            p <- ggplot(df, aes(UMAP_Dim1, UMAP_Dim2)) +
                geom_point(
                    data = df[, 1:2], color = "#e7e7e7",
                    size = point.size
                ) +
                geom_point(
                    size = point.size,
                    aes(color = Clusters)
                    # show.legend = FALSE
                ) +
                facet_wrap(~Clusters) +
                scale_color_manual(
                    values = palette,
                    breaks = gtools::mixedsort(unique(df$Clusters))
                ) +
                labs(
                    x = "UMAP Dim1",
                    y = "UMAP Dim2"
                ) +
                .mytheme +
                theme(
                    legend.text = element_text(size = legend.size)
                ) +
                guides(
                    color = guide_legend(
                        keywidth = 0.4,
                        nrow = 1,
                        override.aes = list(size = 2, stroke = 2)
                    )
                )

            return(p)
        } else {
            df2 <- df %>%
                dplyr::group_by(Clusters) %>%
                dplyr::select(UMAP_Dim1, UMAP_Dim2) %>%
                dplyr::summarize_all(mean)

            p <- ggplot(df,
                aes(UMAP_Dim1, UMAP_Dim2),
                label = TRUE
            ) +
                geom_point(size = point.size, alpha = point.alpha, aes(color = Clusters)) +
                ggrepel::geom_text_repel(
                    data = df2,
                    aes(label = Clusters),
                    color = label.color,
                    bg.color = shadow.color,
                    bg.r = label.shadow,
                    size = label.size,
                    box.padding = 0,
                    seed = 1
                ) +
                scale_color_manual(
                    values = palette,
                    breaks = gtools::mixedsort(df2$Clusters)
                ) +
                guides(
                    color = guide_legend(
                        keywidth = 0.4,
                        override.aes = list(size = 2, stroke = 2)
                    )
                )
        }
    } else {
        mx <- SummarizedExperiment::assay(mtx)
        rownames(mx) <- SummarizedExperiment::rowData(mtx)$name
        df$count <- mx[colorBy, rownames(df)]
        df <- df[order(df$count), ]
        if (log2Norm) df$count <- log2(df$count + 1)
        p <- ggplot(
            df,
            aes(UMAP_Dim1, UMAP_Dim2, label = TRUE)
        ) +
            geom_point(
                size = point.size,
                alpha = point.alpha,
                aes(color = count)
            ) +
            scale_color_gradient2(
                midpoint = max(df$count) * 0.5,
                low = marker.color["low"],
                mid = marker.color["mid"],
                high = marker.color["high"]
            )
        title <- colorBy
    }
    theme <- match.arg(theme)
    if (theme == "mini") {
        p <- p + labs(
            title = title,
            caption = "X-axis: UMAP Dim1\nY-axis: UMAP Dim2"
        ) +
            theme_minimal() +
            theme(
                legend.text = element_text(size = legend.size),
                legend.title = element_blank(),
                axis.title = element_blank(),
                axis.text = element_blank(),
                legend.background = element_rect(
                    color = "#f3eaea"
                ),
                plot.title = element_text(size = title.size)
            )
    } else if (theme == "clear") {
        p <- p + labs(
            title = title,
            x = "UMAP Dim1",
            y = "UMAP Dim2"
        ) +
            tidydr::theme_dr() +
            theme(
                legend.text = element_text(size = legend.size),
                legend.title = element_blank(),
                legend.background = element_rect(
                    color = "#f3eaea"
                ),
                plot.title = element_text(size = title.size),
                panel.grid = element_blank()
            )
    } else {
        stop("theme is required: 'mini' or 'clear'")
    }
    return(p)
}

.mytheme <- theme_minimal() + theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    strip.text = element_blank(),
    panel.spacing = unit(0, "lines"),
    panel.border = element_blank(),
    # panel.spacing.y = unit(-0.5, "lines"),
    panel.grid = element_blank(),
    legend.title = element_blank(),
    legend.position = "top",
    legend.background = element_rect(
        color = "#f3eaea"
    )
)
