library(ArchR)
library(dplyr)
library(fs)
library(networkD3)
addArchRThreads(threads = 32)

basedir <- "/data4/Data_temp/MX/scATAC_220211/MX_test/Result/"
dir_ls(basedir)

pro2 <- loadArchRProject(paste0(basedir, "MX_test2-save/"))
pro4 <- loadArchRProject(paste0(basedir, "MX_test4-save/"))
pro8 <- loadArchRProject(paste0(basedir, "MX_test8-save/"))


proj <- pro8
proj <- addIterativeLSI(
  ArchRProj = proj,
  useMatrix = "TileMatrix",
  name = "IterativeLSI_2",
  iterations = 2,
  clusterParams = list(
    resolution = c(0.8),
    sampleCells = 10000,
    n.start = 10
  ),
  force = TRUE
)
proj <- addClusters(
  input = proj, 
  reducedDims = "IterativeLSI",
  name = "Cluster_iter2",
  force = TRUE
)
proj <- addUMAP(
  ArchRProj = proj,
  reducedDims = "IterativeLSI",
  nNeighbors = 40,
  minDist = 0.4,
  metric = "cosine",
  force = TRUE,
  name = "iter2"
)



T2 <- data.frame(cellNames = paste0(pro2$cellNames), T2 = paste0("T2_", pro2$Clusters), stringsAsFactors = F)
T4 <- data.frame(cellNames = paste0(pro4$cellNames), T4 = paste0("T4_", pro4$Clusters), stringsAsFactors = F)
T8 <- data.frame(cellNames = paste0(pro8$cellNames), T8 = paste0("T8_", pro8$Clusters), stringsAsFactors = F)

tmp <- inner_join(T2, T4)
tmp <- tmp[order(tmp$T4), ]

a <- as.data.frame(table(tmp$T4))[, 2]
tmp$value <- rep(a, a) / 100

links <- data.frame(
  source = tmp$T2,
  target = tmp$T4,
  value = tmp$value
)



nodes <- data.frame(
  name = c(
    as.character(links$source),
    as.character(links$target)
  ) %>% unique()
)


links$IDsource <- match(links$source, nodes$name) - 1
links$IDtarget <- match(links$target, nodes$name) - 1

# Make the Network
p <- sankeyNetwork(
  Links = links, 
  Nodes = nodes,
  Source = "IDsource", 
  Target = "IDtarget",
  Value = "value", 
  NodeID = "name", 
  fontSize = 20,
  sinksRight = FALSE
)


library(htmlwidgets)
saveWidget(p, file = "iteration_24.html")


df <- pro4@embeddings$UMAP$df
colnames(df) <- c("UMAP_Dim1", "UMAP_Dim2")
df$Clusters <- pro4$Clusters

df2 <- df %>%
  group_by(Clusters) %>%
  select(UMAP_Dim1, UMAP_Dim2) %>%
  summarize_all(mean)

ggplot(df,
       aes(UMAP_Dim1, UMAP_Dim2),
       label = T
) +
  geom_point(size = 1, alpha = 0.8, aes(color = Clusters)) +
  labs(color = "Clusters") +
  ggrepel::geom_label_repel(data = df2, aes(label = Clusters)) +
  scale_color_discrete(breaks = gtools::mixedsort(df2$Clusters)) +
  ggtitle("UMAP of IterativeLSI colored by Clusters") +
  theme_classic() +
  theme(
    legend.text = element_text(size = 15),
    legend.title = element_blank(),
    plot.title=element_text(size=20)
  )

source("~/script/tips_func/sankeyForArchR.R")

Name1 <- "Cluster_iter2"
Name2 <- "Cluster_iter4"

p <- Clu.sankey(ArchR.proj = proj, 
                iter.time1 = 2,
                iter.time2 = 4,
                iter.Name1 = Name1, 
                iter.Name2 = Name2
                )
# Make the Network
# library(htmlwidgets)
# saveWidget(p, file = "iteration_24.html")


