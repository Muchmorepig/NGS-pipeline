# rm(list = ls())
# load pkgs ----
options(stringsAsFactors = F)
# library(tidyverse)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(ComplexHeatmap)
library(readr)

# Data preparing ----
rawCounts <- read.table("result/AllCounts.txt", header = T)
rownames(rawCounts) <- rawCounts$Geneid
head(rawCounts)
rawCounts <- rawCounts[, -c(1:6)]
colnames(rawCounts) <- gsub("[.]sort.bam", "", colnames(rawCounts))

(colData <- data.frame(
  name = colnames(rawCounts),
  group = as.factor(rep(c("C1d", "C3d", "W1d", "W3d"), each = 3))
))

rawCounts <- rawCounts[, -c(3:4)]
colData <- colData[-c(3, 4), ]

# cor + PCA -----
dir.create("result/plot", recursive = T)
{ # cor ----
  source("~/script/tips_func/calcHTsize_based_Complexheatmap.R")
  mat <- rawCounts[rowSums(rawCounts > 1) > 1, ]
  mat <- log2(edgeR::cpm(mat) + 1)
  cg <- sort(apply(mat, 1, sd)) %>%
    tail(., 2000) %>%
    names()
  n <- t(scale(t(mat[cg, ])))
  n[n > 2] <- 2
  n[n < -2] <- -2

  p1 <- Heatmap(
    n,
    cluster_columns = F,
    show_row_names = F,
    show_column_dend = F,
    name = "Z-Scores",
    column_title = "SD Top 2000"
  )
  p1
  size1 <- Calc_HTSize(p1)
  size1
  corre <- cor(rawCounts)
  p2 <- Heatmap(
    corre,
    name = "Correlation",
    cluster_rows = F,
    cluster_columns = F,
    show_column_names = F,
    row_names_side = "left",
    heatmap_height = unit(14, "mm") * nrow(corre),
    heatmap_width = unit(15, "mm") * ncol(corre),
    rect_gp = gpar(col = "white"),
    heatmap_legend_param = list(
      direction = "horizontal",
      title_position = "lefttop",
      legend_width = unit(15, "mm") * ncol(corre) * 0.64
    ),
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text(sprintf("%.2f", corre[i, j]), x, y, gp = gpar(fontsize = 10))
    }
  )
  p3 <- draw(p2, heatmap_legend_side = "top")
  size2 <- Calc_HTSize(p3)
  size2


  pdf("result/plot/cor.pdf",
    width = max(size1[1], size2[1]),
    height = max(size1[2], size2[2])
  )
  print(c(p1, p3))
  dev.off()
}
{ # PCA----
  source("~/script/tips_func/PCAplot_based_DESeq2.R")

  dds <- DESeqDataSetFromMatrix(
    countData = rawCounts,
    colData = colData,
    design = ~group
  )
  vsd <- vst(dds)
  PCAplot(vsd)
  PCAplot(vsd, pc = c("PC1", "PC2", "PC3"))
  # htmlwidgets::saveWidget(p,"result_all/plot/PCA-3D.html")

  # pb <- assay(vsd) %>% as_tibble() %>%
  #   tidyr::pivot_longer(cols = everything(),
  #                       names_to = "Sample", values_to = "VSD") %>%
  #   ggplot(data = ., aes(x = Sample, y = VSD)) + geom_boxplot() + theme_bw()

  # pcaData <- PCAplot(vsd,returnData = T)
  # screeplot(pcaData$PCA)
  # pcaData$percentVar %>% data.frame(PC = factor(names(.),levels = names(.)),perc = .) %>%
  #   ggplot(data = .,aes(x = PC, y = perc)) + geom_point(size=2) +
  #   geom_line(size=2,color = "pink",aes(x= c(1:12),y = perc)) +
  #   # geom_smooth(method = "loess",se = FALSE,color = "pink",aes(x= c(1:12),y = perc)) +
  #   theme_bw()
}
# one-one res-----
dds <- DESeq(dds)
count_norm <- counts(dds, normalize = T)
write.csv(count_norm, "result/count_norm.csv")

source("~/script/tips_func/makeCompareInfo.R")
comp <- CompInfo(colData = colData, condition = "group")

comp <- comp[c(1,2,5,6), ]


cutoff_padj <- 0.05
for (x in seq.int(nrow(comp))) {
  # Do the DESeq contrast
  Type <- "group"
  treat <- comp[x, 1] %>% as.character()
  control <- comp[x, 2] %>% as.character()
  message("Trate: ", treat, "\tControl: ", control, "\n")

  res <- results(dds, contrast = c(Type, treat, control), alpha = cutoff_padj)
  suppressMessages(
    res_lfc <- lfcShrink(
      dds = dds, res = res,
      contrast = c(Type, treat, control),
      type = "ashr"
    )
  )
  # output res result
  count_norm_res <- count_norm[, c(
    which(colData[, Type] == treat),
    which(colData[, Type] == control)
  )]

  res_lfc_df <- data.frame(as.data.frame(res_lfc),
    geneId = rownames(res_lfc),
    stringsAsFactors = FALSE
  )
  res_lfc_df <- res_lfc_df[res_lfc_df$baseMean != 0, ]
  write_csv(res_lfc_df,
    file = file.path("result", paste0(comp[x, 3], "_lfc_df.csv"))
  )
  write.csv(count_norm_res,
    file = file.path("result", paste0(comp[x, 3], "_norm_res.csv"))
  )
}

# volcano plot-----
rm(list = ls())
dir.create("result/diff", recursive = T)

gene_alias <- read_csv("~/reference/Mpolymorpha/Mpoly_V5_Tair_anno.csv")

filein <- list.files(path = "result", pattern = "lfc_df.csv$", full.names = T)

for (i in c(1:length(filein))) {
  cat(i, "\n")
  DEG <- read_csv(filein[i], col_types = cols()) %>% na.omit()
  n <- gsub("result/", "", filein[i]) %>% gsub("_lfc_df.csv", "", .)
  gn <- paste0("result/plot/", n, "_volcano.png")
  DEG$Change <- as.factor(
    ifelse(DEG$padj < 0.05 & abs(DEG$log2FoldChange) > 1,
      ifelse(DEG$log2FoldChange > 1, "Up", "Down"),
      "Not"
    )
  )

  DEG %>%
    subset(Change != "Not") %>%
    left_join(., gene_alias, by = "geneId") %>%
    write_csv(paste0("result/diff/", n, "_DIFF.csv"))
  this_title <- paste0(
    n, "  Gene-Up: ",
    nrow(DEG[DEG$Change == "Up", ]),
    "  Gene-Down: ",
    nrow(DEG[DEG$Change == "Down", ]),
    "  Fold Cut-off 1"
  )
  g <- ggplot(
    data = DEG,
    aes(
      x = log2FoldChange,
      y = -log10(pvalue), color = Change
    )
  ) +
    geom_point(aes(size = abs(log2FoldChange), alpha = 0.8)) +
    labs(
      x = "Log2FC",
      y = "-Log10(P-value)",
      title = this_title
    ) +
    scale_size(range = c(1, 6)) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    ggtitle(this_title) +
    # scale_color_manual(values=c("#D6604D", "#d2dae2","#4393C3"))+
    theme_bw() +
    theme(
      plot.title = element_text(size = 15, hjust = 0.05, ),
      legend.background = element_blank(),
      legend.key = element_blank()
    ) +
    scale_color_manual(values = c("blue", "black", "red")) +
    scale_y_continuous(expand = c(0, 2)) +
    guides(
      size = guide_legend(title = "|Log2FC|"),
      color = guide_legend(title = "Sig"),
      alpha = "none"
    )
  ggsave(gn, g, width = 8, height = 8)
}



# go -----
rm(list = ls())
dir.create("result/go")
library(clusterProfiler)
library(magrittr)
library(org.At.tair.db)
library(BuenColors)

filein <- list.files(path = "result/diff", pattern = "csv$", full.names = T)

for (i in c(1:length(filein))) {
  DIFF <- read_csv(filein[i], col_types = cols())
  n <- gsub("result/diff/", "", filein[i]) %>% gsub("_DIFF.csv", "", .)
  gn <- paste0("result/go/", n, "_GO.csv")
  pn <- paste0("result/go/", n, "_GO.png")
  up <- subset(DIFF, Change == "Up")$TAIR %>%
    na.omit() %>%
    unique()
  down <- subset(DIFF, Change == "Down")$TAIR %>%
    na.omit() %>%
    unique()

  dif_list <- list(up = up, down = down)

  cc <- try(compareCluster(
    geneClusters = dif_list,
    fun = "enrichGO",
    OrgDb = "org.At.tair.db",
    keyType = "TAIR",
    ont = "BP",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  ), silent = T)
  if ("try-error" %in% class(cc)) {
    cat(n, "\tNo enrichment found in any of gene cluster\n")
  } else {
    cat(n, "\n")
    gores <- cc@compareClusterResult %>% as.data.frame()
    write_csv(gores, gn)
  }
  p <- dotplot(cc, showCategory = 20, font.size = 10, label_format = 35)
  ggsave(pn, p, height = 10, width = 6)
}
