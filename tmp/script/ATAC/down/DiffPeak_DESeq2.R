# Info ---------------------------
## Script name: DiffPeak_DESeq2.R

## experiment_name:
experiment_name <- "MX"

## Purpose of script:
## 1. using DESeq2 to get diffPeak
## 2. XX

## Author: Guandong Shang
## Date Created: 2021-03-16
## Copyright (c) Guandong Shang, 2020
## Email: shangguandong1996@163.com

# Notes ---------------------------
## 1. XX
## 2. XX   
##

# Prepare -----------------------------------------------------------------

# load up the packages
library(DESeq2)

library(dplyr)
library(ggplot2)
library(tidyr)

library(org.At.tair.db)

# Set Options
options(stringsAsFactors = F)
source("script/function/one_one_parwise.R")


# load up the data
load("result/processed_result/MX_peakRawCount.rda")
peakAnno <- readr::read_csv("result/processed_result/MX_mergePeakAnno.csv")

data_type <- gsub("_rep[123]", "", colnames(peakRawCount_MX))

coldata <- data.frame(row.names = colnames(peakRawCount_MX), 
                      type = data_type, 
                      stringsAsFactors = T)

# creat dir
dir.create("result/Diff_Peak_Anno",recursive = T)
dir.create("result/Diff_Peak_GO_result",recursive = T)
dir.create("plot/MAplot", recursive = T)


## set color map
color_map_cor <- as.character(
  BuenColors::jdb_palette("wolfgang_basic",type = "continuous"))
color_map_tissue <- BuenColors::jdb_palette("brewer_spectra")[c(1,2,3,4,6,7,8)]

# set cutoff
cutoff_padj <- 0.05
cutoff_lfc <- 1


# DESeq2 set ------------------------------------------------------------------
# dds set
dds <- DESeqDataSetFromMatrix(countData = peakRawCount_MX, 
                              colData = coldata, 
                              design = ~ type)

# vsd & PCA
vsd <- vst(dds,blind = FALSE)
assign(paste0("vsd_", experiment_name), assay(vsd))
save(list = paste0(paste0("vsd_", experiment_name)), 
     file = paste0("result/processed_result/", experiment_name, "_vsd.rda"))

pcaData <- plotPCA(vsd, intgroup=c("type"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pcaData %>% 
  as_tibble() %>% 
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(color = type), size = 4) +
  theme(legend.key = element_rect(fill = "transparent", colour = NA)) +
  scale_color_manual(values = color_map_tissue) +
  ggrepel::geom_text_repel(aes(label = name)) +
  theme_bw() +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  ggsave(paste0("plot/", experiment_name, "_vsdPCA.pdf"),
         width = 6, height = 4)



# DESeq2 one vs one -------------------------------------------------------

# dds one vs one 
## although we can use the above result through 
## results(dds_LRT, alpha = 0.05, test="Walt", name = ……)
dds <- DESeq(dds)

count_norm <- counts(dds, normalized = TRUE)
assign(paste0("peakCountNorm_", experiment_name), count_norm)
save(list = paste0("peakCountNorm_", experiment_name), 
     file = paste0("result/processed_result/", experiment_name, "_peakNormCount.rda"))

# compare expand
compare_expand <- expand.grid(unique(data_type), 
                              unique(data_type), 
                              stringsAsFactors = F) %>% 
  slice(2, 3) %>% 
  mutate(key = paste0(Var1, "_VS_", Var2))


# get all results
lapply(1:dim(compare_expand)[1], one_one_parwise, "type", dds)

