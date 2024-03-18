if (!"dplyr" %in% (.packages())) suppressMessages(library(dplyr))
if (!"ggplot2" %in% (.packages())) suppressMessages(library(ggplot2))
if (!"gtools" %in% installed.packages()[, 1]) install.packages("gtools")
if (!"paletteer" %in% installed.packages()[, 1]) install.packages("paletteer")


umapForArchR <- function(ArchR.proj = proj,
                         embed.umap = "UMAP_Tile_8",
                         embed.cluster = "Cluster_Tile_8",
                         colors = paletteer::paletteer_d("ggsci::category20c_d3"),
                         label.size = 3.6,
                         label.shadow = 0.25,
                         label.color = "#000103",
                         shadow.color = "#f5f0f0",
                         point.size = 0.8,
                         point.alpha = 0.8,
                         theme = c("mini", "clear"),
                         title = NULL,
                         title.size = 13,
                         legend.size = 10) {
    df <- ArchR.proj@embeddings[[embed.umap]]$df
    colnames(df) <- c("UMAP_Dim1", "UMAP_Dim2")

    if (stringr::str_to_title(tolower(embed.cluster)) == "Sample") {
        df$Clusters <- as.vector(proj@cellColData[, embed.cluster])
    } else {
        df$Clusters <- proj@cellColData[, embed.cluster]
    }

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
            values = colors,
            breaks = gtools::mixedsort(df2$Clusters)
        )
    # scale_color_discrete()

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
            ) +
            guides(
                color = guide_legend(
                    keywidth = 0.4,
                    override.aes = list(size = 2, stroke = 2)
                )
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
            ) +
            guides(
                color = guide_legend(
                    keywidth = 0.4,
                    override.aes = list(size = 2, stroke = 2)
                )
            )
    } else {
        stop("theme must be 'mini' or 'clear'")
    }
    p
}
