library(ggplot2)
library(ggraph)
library(igraph)

gobarplot <- function(dat = g, size = 8, tt = "") {
    if (!"BuenColors" %in% (.packages())) library(BuenColors)
    # dat <- dplyr::arrange(dat, -log10(pvalue))
    dat <- dplyr::arrange(dat, -log10(p.adjust))
    dat$Description <- stringr::str_to_title(
        factor(dat$Description, levels = dat$Description)
    )
    ggplot(
        dat,
        aes(
            x = Description,
            # y = -log10(pvalue),
            y = -log10(p.adjust),
            label = Description
        )
    ) +
        geom_bar(
            stat = "identity",
            position = position_stack(reverse = T),
            fill = jdb_color_map(c("CLP")), alpha = 0.75
        ) +
        theme_minimal() +
        coord_flip() +
        L_border() +
        # labs(x = "Pathway Name", y = "-log10(pvalue)", title = paste0("GO ENRICH: ", tt)) +
        labs(x = "Pathway Name", y = "-log10(padjust)", title = paste0("GO ENRICH: ", tt)) +
        theme(legend.position = "none", axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
        geom_text(color = "black", position = position_stack(vjust = 0.02), hjust = 0, size = size) +
        scale_y_continuous(expand = c(0, 0.01)) +
        theme(
            title = element_text(size = 15),
            text = element_text(size = 15),
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 15)
        )
}


gonetplort <- function(x,
                       showCategory = 5,
                       # foldChange = NULL,
                       layout = "kk",
                       alias = read.table(
                           "/data3/wanglab/wmc/reference/Oryza_sativa/annotation/OS_symbol_simple2.txt",
                           header = TRUE
                       ),
                       colorEdge = FALSE,
                       circular = FALSE,
                       node_label = "all",
                       cex_category = 1,
                       cex_gene = 2,
                       cex_label_category = 1,
                       cex_label_gene = 1,
                       color_category = "#E5C494",
                       color_gene = "#B3B3B3",
                       ...) {
    label_size_category <- 5
    label_size_gene <- 5

    node_label <- match.arg(node_label, c("category", "gene", "all", "none"))
    alias <- alias
    if (circular) {
        layout <- "linear"
        geom_edge <- ggraph::geom_edge_arc
    } else {
        geom_edge <- ggraph::geom_edge_link
    }

    shadowtext_category <- shadowtext_gene <- TRUE

    geneSets <- .extract_geneSets(x, showCategory, alias)

    g <- .list2graph(geneSets)

    # if (!inherits(x, "list")) {
    #     foldChange <- .fc_readable(
    #         x # , foldChange
    #     )
    # }

    size <- sapply(geneSets, length)

    igraph::V(g)$size <- min(size) / 2
    n <- length(geneSets)
    igraph::V(g)$size[1:n] <- size
    if (colorEdge) {
        igraph::E(g)$category <- rep(names(geneSets), sapply(geneSets, length))
        edge_layer <- geom_edge(ggplot2::aes_(color = ~category), alpha = .8)
    } else {
        edge_layer <- geom_edge(alpha = .8, colour = "darkgrey")
    }

    # if (!is.null(foldChange)) {
    #     fc <- foldChange[V(g)$name[(n + 1):length(V(g))]]
    #     V(g)$color <- NA
    #     V(g)$color[(n + 1):length(V(g))] <- fc
    #     show_legend <- c(TRUE, FALSE)
    #     names(show_legend) <- c("color", "size")
    #     p <- ggraph::ggraph(g, layout = layout, circular = circular)
    #     p <- p + edge_layer +
    #         geom_node_point(aes_(color = ~ I(color_category), size = ~size),
    #             data = p$data[1:n, ]
    #         ) +
    #         scale_size(range = c(3, 8) * cex_category) +
    #         ggnewscale::new_scale_color() +
    #         geom_node_point(aes_(color = ~ as.numeric(as.character(color)), size = ~ I(3 * cex_gene)),
    #             data = p$data[- (1:n), ], show.legend = show_legend
    #         ) +
    #         scale_colour_gradient2(
    #             name = "fold change", low = "blue",
    #             mid = "white", high = "red",
    #             guide = guide_colorbar(order = 2)
    #         )
    # } else {
    #     igraph::V(g)$color <- color_gene
    #     igraph::V(g)$color[1:n] <- color_category
    #     p <- ggraph::ggraph(g, layout = layout, circular = circular)
    #     p <- p + edge_layer +
    #         geom_node_point(aes_(color = ~ I(color), size = ~size), data = p$data[1:n, ]) +
    #         scale_size(range = c(3, 8) * cex_category) +
    #         geom_node_point(aes_(color = ~ I(color), size = ~ I(3 * cex_gene)),
    #             data = p$data[- (1:n), ], show.legend = FALSE
    #         )
    # }
    igraph::V(g)$color <- color_gene
    igraph::V(g)$color[1:n] <- color_category
    p <- ggraph::ggraph(g, layout = layout, circular = circular)
    p <- p + edge_layer +
        geom_node_point(aes_(color = ~ I(color), size = ~size), data = p$data[1:n, ]) +
        scale_size(range = c(3, 8) * cex_category) +
        geom_node_point(aes_(color = ~ I(color), size = ~ I(3 * cex_gene)),
            data = p$data[-(1:n), ], show.legend = FALSE
        )
    p <- p + ggplot2::theme_void()
    if (node_label == "category") {
        p <- .add_node_label(
            p = p, data = p$data[1:n, ], label_size_node = label_size_category,
            cex_label_node = cex_label_category, shadowtext = shadowtext_category
        )
    } else if (node_label == "gene") {
        p <- .add_node_label(
            p = p, data = p$data[-c(1:n), ], label_size_node = label_size_gene,
            cex_label_node = cex_label_gene, shadowtext = shadowtext_gene
        )
    } else if (node_label == "all") {
        p <- .add_node_label(
            p = p, data = p$data[-c(1:n), ], label_size_node = label_size_gene,
            cex_label_node = cex_label_gene, shadowtext = shadowtext_gene
        )
        p <- .add_node_label(
            p = p, data = p$data[1:n, ], label_size_node = label_size_category,
            cex_label_node = cex_label_category, shadowtext = shadowtext_category
        )
    }
    # if (!is.null(foldChange)) {
    #     p <- p + guides(
    #         size = guide_legend(order = 1),
    #         color = guide_colorbar(order = 2)
    #     )
    # }
    return(p)
}




.geneInCategory <- function(x, alias) {
    if (is.null(alias)) {
        setNames(strsplit(x$geneID, "/", fixed = TRUE), rownames(x))
    } else {
        setNames(
            lapply(
                strsplit(x$geneID, "/", fixed = TRUE),
                function(i) {
                    alias[which(alias$RAP_ID %in% i), 2]
                }
            ),
            x$ID
        )
    }
    # setNames(strsplit(x$geneID, "/", fixed = TRUE), rownames(x))
}

.update_n <- function(x, showCategory) {
    if (!is.numeric(showCategory)) {
        if (inherits(x, "list")) {
            showCategory <- showCategory[showCategory %in% names(x)]
        }
        return(showCategory)
    }
    n <- showCategory
    if (inherits(x, "list")) {
        nn <- length(x)
    } else {
        nn <- nrow(x)
    }
    if (nn < n) {
        n <- nn
    }
    return(n)
}

.list2graph <- function(inputList) {
    x <- .list2df(inputList)
    g <- igraph::graph.data.frame(x, directed = FALSE)
    return(g)
}

.extract_geneSets <- function(x, n, alias) {
    n <- .update_n(x, n)

    if (inherits(x, "list")) {
        geneSets <- x
    } else {
        geneSets <- .geneInCategory(x, alias)
        # geneSets <- DOSE::geneInCategory(x) ## use core gene for gsea result
        y <- as.data.frame(x)
        geneSets <- geneSets[y$ID]
        names(geneSets) <- stringr::str_to_title(y$Description)
    }
    if (is.numeric(n)) {
        return(geneSets[1:n])
    }
    return(geneSets[n]) ## if n is a vector of Description
}

# Convert a list of gene IDs to data.frame object.
.list2df <- function(inputList) {
    # ldf <- lapply(1:length(inputList), function(i) {
    ldf <- lapply(seq_len(length(inputList)), function(i) {
        data.frame(
            categoryID = rep(
                names(inputList[i]),
                length(inputList[[i]])
            ),
            Gene = inputList[[i]]
        )
    })

    do.call("rbind", ldf)
}


# .fc_readable <- function(x, foldChange = NULL) {
#     if (is.null(foldChange)) {
#         return(NULL)
#     }

#     if (x@readable) {
#         gid <- names(foldChange)
#         if (is(x, "gseaResult")) {
#             ii <- gid %in% names(x@geneList)
#         } else {
#             ii <- gid %in% x@gene
#         }
#         gid[ii] <- x@gene2Symbol[gid[ii]]
#         names(foldChange) <- gid
#     }
#     return(foldChange)
# }


.add_node_label <- function(p, data, label_size_node, cex_label_node, shadowtext) {
    segment.size <- .get_ggrepel_segsize()
    if (shadowtext) {
        p <- p + geom_node_text(aes_(label = ~name),
            data = data,
            size = label_size_node * cex_label_node, bg.color = "white",
            repel = TRUE, segment.size = segment.size
        )
    } else {
        p <- p + geom_node_text(aes_(label = ~name),
            data = data,
            size = label_size_node * cex_label_node, repel = TRUE,
            segment.size = segment.size
        )
    }
    return(p)
}

.get_ggrepel_segsize <- function(default = 0.2) {
    getOption("ggrepel.segment.size", default = default)
}



# compClusterResult <- function(x) {
#     x <- as.data.frame(x)
#     reslist <- split(x, x$Cluster)
#     if ("core_enrichment" %in% colnames(x)) {
#         res <- lapply(reslist, function(y) {
#             setNames(
#                 strsplit(y$core_enrichment, split = "/", fixed = TRUE),
#                 y$ID
#             )
#         })
#     } else {
#         res <- lapply(reslist, function(y) {
#             setNames(
#                 strsplit(y$geneID, split = "/", fixed = TRUE),
#                 y$ID
#             )
#         })
#     }
#     res[vapply(res, length, numeric(1)) != 0]
# }