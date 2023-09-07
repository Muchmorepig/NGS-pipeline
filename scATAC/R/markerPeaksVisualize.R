library(dplyr)
library(readxl)
library(ggplot2)
library(data.table)
library(patchwork)
library(rtracklayer)
library(aplot)

load("dataSheet/MSU7.Rdata")
geneAnnotation$genes
geneAnnotation$exons



marker_file <- "./Result/Merge/Marker/MarkerPeaks.xlsx"

sheets <- excel_sheets(marker_file)
markerList <- purrr::map(
    sheets, ~ read_excel(marker_file, sheet = .)
)
names(markerList) <- sheets
sapply(markerList, nrow)





locus <- fread("dataSheet/IRGSP_OsLocus.csv")

(gene <- locus[ID == "Os02g0747400"])
unique(gene$type)

sel <- gene[type == "gene"]
pmyc <- ggplot(sel, aes(group = ID)) +
    geom_segment(
        data = sel,
        aes(
            x = min(start), y = 0.02,
            xend = max(end) + 50, yend = 0.02
        ),
        arrow = arrow(length = unit(0.15, "cm"), type = "closed")
    ) +
    geom_rect(
        aes(xmin = start, xmax = end, ymin = 0, ymax = 0.04),
        fill = "#6db8f6"
    ) +
    theme_void() +
    theme(aspect.ratio = 0.05) +
    geom_text(
        data = data.frame(
            x = sel$start + (sel$end - sel$start) / 2,
            y = 0.02,
            ID = "Os02g0747400",
            text = "Os02g0747400"
        ),
        aes(x = x, y = y, label = text),
        size = 2.5
    )



bw <- import.bw("Result/Merge/GroupBigWigs/Cluster_Peaks/C1-TileSize-40-normMethod-ReadsInTSS-ArchR.bw") %>%
    as.data.table() %>%
    select(-width, -strand)

colnames(bw)

long_bw <- melt(bw, id.vars = c("seqnames", "score")) %>%
    select(-variable) %>%
    arrange(seqnames, value)



tmp <- markerList$C1 %>%
    select(-FDR, -MeanDiff)
range <- subset(tmp, RAP_ID == "Os02g0747400")

plot_range <- long_bw %>%
    filter(seqnames == unique(range$seqnames) &
        value <= max(range$end) &
        value >= min(range$start))

ptest <- ggplot(plot_range, aes(x = value, y = score)) +
    geom_area(fill = "#f26868") +
    theme_bw(base_size = 16) +
    scale_y_continuous(
        expand = c(0, 0),
        n.breaks = 6,
        limits = c(0, max(plot_range$score * 1.05))
    ) +
    theme(
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid = element_blank(),
        axis.ticks.length.y = unit(0.25, "cm")
    )

ptest %>% insert_bottom(pmyc)


file_loc <- list.files("Result/Merge/GroupBigWigs/Cluster_Peaks",
    pattern = "bw$",
    full.names = TRUE
)

library(fs)
combres <- lapply(file_loc, function(file) {
    print(file)
    bw <- import.bw(file) %>%
        as.data.table() %>%
        select(-width, -strand)
    bw$sample_name <- gsub(
        pattern = "\\-[0-9]*|\\-[A-z]*|.bw",
        replacement = "",
        x = path_file(file)
    )
    return(bw)
}) %>%
    do.call(rbind, .) %>%
    as.data.table()


table(combres$sample_name)

st_time <- Sys.time()
long_bw_com <- melt(
    combres,
    id.vars = c("seqnames", "score", "sample_name")
) %>%
    select(-variable) %>%
    arrange(sample_name, seqnames, value)

ed_time <- Sys.time()
ed_time - st_time

long_bw_com$seqnames <- unfactor(long_bw_com$seqnames)

plot_range <- long_bw_com %>%
    filter(seqnames == unique(range$seqnames) &
        value <= max(range$end) &
        value >= min(range$start))


pcomb <- ggplot(plot_range, aes(x = value, y = score)) +
    geom_area(fill = "#f26868", show.legend = FALSE) +
    theme_bw() +
    scale_y_continuous(
        expand = c(0, 0),
        n.breaks = 5,
        limits = c(0, max(plot_range$score * 1.05))
    ) +
    # scale_fill_brewer(palette = "Set1", name = "")
    theme(
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid = element_blank(),
        strip.placement = "outside",
        strip.background = element_rect(color = "white"),
        panel.spacing = unit(0, "mm"),
        axis.ticks.length.y = unit(2, "mm")
    ) +
    facet_wrap(~sample_name, ncol = 1, strip.position = "left") +
    ylab("Normalized densitys")

p1 <- pcomb %>% insert_bottom(pmyc, height = 0.04)


library(ArchR)

proj <- loadArchRProject("Result/Merge")

p2 <- plotEmbedding(
    ArchRProj = proj,
    colorBy = "GeneScoreMatrix",
    name = "Os02g0747400",
    embedding = "UMAP_Peaks"
)

pname <- paste0(
    getOutputDirectory(proj),
    "/Plots/markerPlot/",
    "Os02g0747400.pdf"
)

pdf(pname, width = 4, height = 10)
p1
dev.off()
