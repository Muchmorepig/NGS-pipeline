if (!"ComplexHeatmap" %in% (.packages())) suppressMessages(library(ComplexHeatmap))
if (!"magrittr" %in% (.packages())) library(magrittr)
dotplotForArchR <- function(geneScoreMatrix = GSmat,
                            name = "Cluster_Tile_6",
                            genes = NULL,
                            showclusters = NULL,
                            col.min = -2.5,
                            col.max = 2.5,
                            scale = TRUE,
                            threshold = 0,
                            dotsize = 10,
                            legend.dotsize = dotsize - 1,
                            geneFontsize = 8,
                            dotcolor = c("#2f0154", "#19b99f", "#f58114"),
                            rowsp = NULL,
                            borderlwd = 2.5,
                            ...) {
    if (!methods::is(geneScoreMatrix, "SummarizedExperiment")) {
        stop("需要 geneScoreMatrix 是 SummarizedExperiment object")
    }

    if (!is.null(rowsp)) {
        rowsp <- rowsp[!duplicated(genes)]
    }

    genes <- genes[!duplicated(genes)]

    mat <- SummarizedExperiment::assay(geneScoreMatrix)
    rownames(mat) <- SummarizedExperiment::rowData(geneScoreMatrix)$name
    submat <- mat[genes, ]

    clus <- data.frame(
        "Clusters" = geneScoreMatrix@colData[name][, 1],
        "cellNames" = rownames(geneScoreMatrix@colData[name])
    )

    if (!is.null(showclusters)) {
        clus <- clus[which(clus$Clusters %in% showclusters), ]
    }

    plotdata <- .plotData(
        submat, clus, genes, col.min,
        col.max, scale, threshold
    )
    perc_mat <- plotdata$percent
    exp_scaled <- plotdata$scaled


    cluster_anno <- colnames(exp_scaled)
    # anno_color <- .col[sample(1:28)[seq_along(cluster_anno)]]
    anno_color <- .selcolor(cluster_anno)
    if (length(unique(clus$Clusters)) < 8) {
        ncol <- 1
    } else {
        ncol <- 2
    }

    column_ha <- HeatmapAnnotation(
        Cluster = cluster_anno,
        col = list(
            Cluster = setNames(anno_color, unique(cluster_anno))
        ),
        na_col = "grey",
        show_annotation_name = F,
        annotation_legend_param = list(
            Cluster = list(
                # direction = "horizontal",
                ncol = ncol
            )
        )
    )

    col_fun <- circlize::colorRamp2(
        c(min(exp_scaled), 0, max(exp_scaled)),
        dotcolor
        # viridis(20)[c(1, 10, 20)]
    )

    layer_fun <- function(j, i, x, y, w, h, fill) {
        grid.rect(
            x = x, y = y, width = w, height = h,
            gp = gpar(col = NA, fill = NA)
        )
        grid.points(
            x = x, y = y,
            gp = gpar(col = col_fun(pindex(exp_scaled, i, j))),
            size = pindex(perc_mat, i, j) * unit(dotsize, "mm"),
            pch = 19
        )
    }

    lgd_list <- list(
        Legend(
            col_fun = col_fun,
            title = "Average\nAccessibility",
            grid_height = unit(40, "mm"),
            grid_width = unit(8, "mm")
            # direction = "horizontal",
            # title_position = "topcenter",
            # legend_height = unit(50, "mm")
        ),
        Legend(
            labels = c(20, 40, 60, 80),
            title = "Accessible%",
            type = "points",
            background = NULL,
            pch = 19,
            # nrow = 1,
            size = c(0.20, 0.40, 0.60, 0.80) * unit(legend.dotsize, "mm"),
            grid_height = unit(legend.dotsize, "mm"),
            grid_width = unit(legend.dotsize, "mm")
            # direction = "horizontal",
            # title_position = "topcenter",
        )
    )

    if (!is.null(rowsp)) {
        rowsp <- factor(rowsp, unique(rowsp))
    }

    hp <- Heatmap(
        exp_scaled,
        column_title = "",
        col = col_fun,
        layer_fun = layer_fun,
        top_annotation = column_ha,
        rect_gp = gpar(type = "none"),
        row_names_gp = gpar(fontsize = geneFontsize),
        border = TRUE,
        border_gp = gpar(col = "black", lwd = borderlwd),
        column_names_rot = 0,
        column_names_side = "top",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        show_heatmap_legend = FALSE,
        column_names_centered = TRUE,
        row_split = rowsp,
        row_gap = unit(0, "mm")
    )

    ComplexHeatmap::draw(
        hp,
        annotation_legend_list = lgd_list,
        annotation_legend_side = "left"
    )
}


.MinMax <- function(data, min, max) {
    data[data > max] <- max
    data[data < min] <- min
    return(data)
}

.PercentAbove <- function(x, threshold) {
    return(length(x[x > threshold]) / length(x))
}

.plotData <- function(submat, clus, genes, col.min, col.max, scale, threshold) {
    cc <- unique(clus$Clusters)
    cc <- cc[gtools::mixedorder(cc)]
    perc_mat <- lapply(cc, function(i) {
        cn <- clus[clus$Clusters == i, "cellNames"]
        dat.use <- submat[, cn]
        per <- apply(dat.use, 1, .PercentAbove, threshold = 0)
        per <- data.frame(per)
        colnames(per) <- i
        return(per)
    }) %>%
        do.call(cbind, .) %>%
        as.matrix()


    exp_mat <- lapply(cc, function(i) {
        cn <- clus[clus$Clusters == i, "cellNames"]
        dat.use <- submat[, cn]
        avgexp <- apply(dat.use, 1, function(x) {
            return(mean(expm1(x)))
        })
        avgexp <- data.frame(avgexp, row.names = names(avgexp))
        colnames(avgexp) <- i
        return(avgexp)
    }) %>% do.call(cbind, .)

    exp_scaled <- sapply(
        genes,
        function(x) {
            dat.use <- as.numeric(exp_mat[x, ])
            if (scale) {
                dat.use <- scale(x = dat.use)
                dat.use <- .MinMax(data = dat.use, min = col.min, max = col.max)
            } else {
                dat.use <- log1p(x = dat.use)
            }
            return(dat.use)
        }
    ) %>% t()
    colnames(exp_scaled) <- colnames(exp_mat)

    return(list("percent" = perc_mat, "scaled" = exp_scaled))
}


.selcolor <- function(cluster_anno) {
    col <- c(
        "#cc1616", "#d0595a", "#f786a8", "#921f15",
        "#f7b090", "#ff7f44", "#f85408", "#f01778",
        "#b7889d", "#ecd0d0", "#17e7ee", "#8dcaec",
        "#04a5fc", "#66c2a5", "#3bd608", "#a9e097",
        "#2c6917", "#cab9f5", "#df8633", "#05968e",
        "#5e6166", "#1872a3", "#7774c2", "#393b78",
        "#0b13f1", "#a00b98", "#63065d", "#2e012e"
    )
    return(sample(col)[seq_along(cluster_anno)])
}