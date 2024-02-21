# TODO: 统计每个cluster的样本占比，绘图
suppressMessages(library(ArchR))
library(tidyr)
library(dplyr)
library(ggplot2)

# 加载项目-----------
proj <- "Result/Merge"
proj <- loadArchRProject(proj)

# 统计cluster数目和样本占比--------------
cellcol <- proj@cellColData %>%
    as_tibble() %>%
    select(Sample, Cluster_Peaks) %>%
    nest(data = Sample) %>%
    rowwise(data) %>%
    mutate(
        cell_count = nrow(data),
        SIM16d = 100 * nrow(data[data == "SIM16d", ]) / nrow(data),
        SIM12d = 100 * nrow(data[data == "SIM12d", ]) / nrow(data),
        CIM22d = 100 * nrow(data[data == "CIM22d", ]) / nrow(data)
    )

# 绘图 ---------------
# Cluster细胞数目
p1 <- cellcol %>%
    select(Cluster_Peaks, cell_count) %>%
    ggplot(aes(x = Cluster_Peaks, y = cell_count)) +
    geom_bar(stat = "identity") +
    scale_y_continuous(expand = c(0, 10)) +
    labs(x = "Clusters", y = "Cell Count") +
    ggprism::theme_prism() +
    theme(
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank()
    ) +
    coord_flip() +
    geom_text(
        data = data.frame(
            Cluster_Peaks = "C1",
            cell_count = 1400,
            lab = paste0("Total: ", sum(cellcol$cell_count))
        ),
        aes(label = lab),
        size = 6
    )

# 每个Cluster的样本占比
pdat <- tibble(
    cluster_peaks = cellcol$Cluster_Peaks,
    CIM22d = cellcol$CIM22d,
    SIM12d = cellcol$SIM12d,
    SIM16d = cellcol$SIM16d
) %>%
    pivot_longer(
        cols = CIM22d:SIM16d,
        names_to = "Period", values_to = "Proportion"
    ) %>%
    mutate(pos = Proportion * 0.5)

p2 <- ggplot(pdat, aes(x = cluster_peaks, y = Proportion, fill = Period)) +
    geom_bar(position = "stack", stat = "identity") +
    scale_y_continuous(expand = c(0, 1)) +
    scale_fill_manual(values = c("#fa7d54", "#f8a897", "#fbd8cf")) +
    labs(x = "Clusters", y = "Proportion") +
    ggprism::theme_prism() +
    theme(
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank()
    ) +
    coord_flip()


# CIM22d->SIM12d->SIM16d, 样本 cluster 比例变化

coldatBySamples <- proj@cellColData %>%
    as_tibble() %>%
    select(Sample, Cluster_Peaks) %>%
    nest(cp = Cluster_Peaks) %>%
    rowwise() %>%
    mutate(pp = list(as_tibble(
        prop.table(table(cp)) * 100,
        n = "proportion"
    )))

coldatBySamples

props <- do.call(rbind, coldatBySamples$pp) %>%
    mutate(
        cc = rep(1:3, sapply(coldatBySamples$pp, nrow))
    )

View(props)


p3 <- ggplot(props, aes(x = cc, y = proportion, fill = Cluster_Peaks)) +
    geom_area() +
    scale_fill_manual(
        values = paletteer::paletteer_d("awtools::bpalette"),
        breaks = gtools::mixedsort(props$Cluster_Peaks)
    ) +
    scale_y_continuous(
        expand = c(0, 0.05),
        breaks = c(25, 50, 75),
        labels = c("25%", "50%", "75%")
    ) +
    scale_x_continuous(
        expand = c(0, 0),
        breaks = c(1, 2, 3),
        labels = coldatBySamples$Sample
    ) +
    ggprism::theme_prism() +
    theme(
        axis.title = element_blank(),
        # axis.text.y = element_blank(),
        # axis.ticks.y = element_blank(),
        axis.line.y = element_blank()
    )


# 拼图, 保存---------
library(patchwork)

pdf(paste0(getOutputDirectory(proj), "/Plots/Merge_cluster_proportion.pdf"), width = 21)
p1 + p2 + p3
dev.off()
