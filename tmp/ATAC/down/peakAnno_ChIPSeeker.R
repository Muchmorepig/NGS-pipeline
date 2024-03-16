# Info ---------------------------
## Script name: peakAnno_ChIPSeeker.R

## experiment_name:
experiment_name <- "MX"

## Purpose of script:
## 1. using ChIPSeeker to annotate Peak
## 2. XX

## Author: Guandong Shang
## Date Created: 2021-04-02
## Copyright (c) Guandong Shang, 2021
## Email: shangguandong1996@163.com

# Notes ---------------------------
## 1. XX
## 2. XX   
##

# Prepare -----------------------------------------------------------------

# load up the packages
library(ChIPseeker)


# Set Options
options(stringsAsFactors = F)

# load up the data
Txdb <- AnnotationDbi::loadDb("/data/sgd_data/reference/annoation/Athaliana/Araport11/Txdb_gff_Araport11.sqlite")
load("result/processed_result/MX_mergePeakGR.rda")

gene_alias <- readr::read_csv("/data/sgd_data/reference/annoation/Athaliana/Araport11/gene_aliases_20190630_tidy_full_description.csv")


# peak Anno ---------------------------------------------------------------

peakAnno <- annotatePeak(peak = mergePeak_GR_MX,
                         tssRegion = c(-500, 500),
                         TxDb = Txdb,
                         level = "gene")


# outputPeak --------------------------------------------------------------

library(dplyr)

peakAnno %>% 
  as_tibble() %>% 
  mutate(geneChr = paste0("Chr", geneChr)) %>% 
  left_join(gene_alias, by = c("geneId" = "name")) %>% 
  readr::write_csv(file = paste0("result/processed_result/", experiment_name, "_mergePeakAnno.csv"))
