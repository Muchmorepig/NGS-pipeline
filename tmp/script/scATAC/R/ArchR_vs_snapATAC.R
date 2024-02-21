# Comparsion with snapATAC ------------------------------------------------
df <- merge(df1, df2, by = "barcode")
snap_cluster <- read.table("./snap_cluster.txt", header = TRUE)

df <- merge(df, snap_cluster, by = "barcode")

table(df$cluster.x, df$cluster.y)

cM <- confusionMatrix(paste0(df$cluster.x), paste0(df$cluster.y))
cM <- cM / Matrix::rowSums(cM)
p <- pheatmap::pheatmap(
    mat = as.matrix(cM),
    color = paletteContinuous("whiteBlue"),
    border_color = "black"
)
p


cM <- confusionMatrix(paste0(df$cluster.x), paste0(df$cluster))
cM <- cM / Matrix::rowSums(cM)
p <- pheatmap::pheatmap(
    mat = as.matrix(cM),
    color = paletteContinuous("whiteBlue"),
    border_color = "black"
)
p


plot_df <- proj2@embeddings$UMAP$df
colnames(plot_df) <- c("UMAP1", "UMAP2")
plot_df$barcode <- substring(row.names(plot_df), 7)


plot_df <- merge(plot_df, df, by = "barcode")
plot_df$cluster <- paste0("C-", plot_df$cluster)

p.x <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2)) +
    geom_point(aes(color = cluster.x), size = 0.7) +
    scale_color_manual(values = BuenColors::jdb_palette("corona")) +
    theme_classic()
p.x

p.y <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2)) +
    geom_point(aes(color = cluster.y), size = 0.7) +
    scale_color_manual(values = BuenColors::jdb_palette("corona")) +
    theme_classic()
p.y



p.z <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2)) +
    geom_point(aes(color = cluster), size = 0.7) +
    # scale_color_manual(values = BuenColors::jdb_palette("corona")) +
    theme_classic()
p.z








snap_cluster <- read.table("./snap_cluster.txt", header = TRUE)
plot_df <- merge(plot_df, snap_cluster, by.x = "barcode", by.y = "barcode")
plot_df$cluster <- factor(plot_df$cluster)

plot_df2 <- plot_df

table(plot_df$color, plot_df$cluster)



pp <- ggplot(plot_df2, aes(x = x, y = y)) +
    geom_point(aes(color = cluster), size = 0.7) +
    # scale_color_manual(values = BuenColors::jdb_palette("corona")) +
    theme_classic()
pp