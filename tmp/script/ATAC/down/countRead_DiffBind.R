# Info ---------------------------
## Script name: countRead_DiffBind.R

## experiment_name:
experiment_name <- ""

## Purpose of script:
## 1. using DiffBind to mergePeak and countReads
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
library(DiffBind)

# Set Options
options(stringsAsFactors = F)
dir.create("result/processed_result", recursive = T)

# load up the data
bam_files <- list.files("rawdata/bam", pattern = "bam$", full.names = TRUE)
peak_files <- list.files("rawdata/peak", full.names = TRUE)

tissue <- rep(c("Col_0", "SPY", "UBR7"), each = 2)

sample_info <- data.frame(
  SampleID = paste(tissue, c("rep1", "rep2"), sep = "_"),
  
  Tissue = tissue,
  
  Replicate = 1:2,
  
  bamReads = bam_files, 
  
  Peaks = peak_files,
  
  PeakCaller = "narrow",
  
  stringsAsFactors = F
  
)  

## set color map
color_map_cor <- as.character(
  BuenColors::jdb_palette("wolfgang_basic",type = "continuous"))
color_map_tissue <- BuenColors::jdb_palette("brewer_spectra")[c(1,2,3,4,6,7,8)]

# Step1 importing data ----------------------------------------------------
dba_meta <- dba(minOverlap = 2, sampleSheet = sample_info)
dba_meta


# Step2 counting reads ----------------------------------------------------
dba_count <- dba.count(dba_meta, minOverlap = 2)
dba_count

assign(paste0("dba_count_", experiment_name), dba_count)
save(list = paste0("dba_count_", experiment_name), 
     file = paste0("result/processed_result/", experiment_name, "_dbaCount.rda"))


pdf(paste0("plot/",experiment_name,"_Sample_Cor.pdf"),
    width = 8,height = 7)
# 样本相关性
dba.plotHeatmap(dba_count,
                colScheme = color_map_cor,
                RowAttributes = DBA_TISSUE,
                ColAttributes = F,
                rowSideCols = color_map_tissue)

dba.plotPCA(dba_count,
            vColors = color_map_tissue,
            label = DBA_ID)

# peak之间的相关性,可以用来看看组织特异的peak
dba.plotHeatmap(dba_count,
                colScheme = color_map_cor,
                score = DBA_SCORE_TMM_READS_FULL,
                ColAttributes = DBA_TISSUE,
                colSideCols = color_map_tissue,
                correlations = F)
dev.off()

# merge_peak & rawcount ---------------------------------------------------

# Rawcount
peak_count_list <- dba_count$peaks
names(peak_count_list) <- dba_count$samples$SampleID

Peak_count <- lapply(peak_count_list, function(x) {x$Reads})
ATAC_Peak_count <- do.call(cbind, Peak_count)
rownames(ATAC_Peak_count) <- paste0(experiment_name,"_",1:dim(ATAC_Peak_count)[1])

assign(paste0("peakRawCount_", experiment_name), ATAC_Peak_count)
save(list = paste0("peakRawCount_", experiment_name),
     file = paste0("result/processed_result/", experiment_name, "_peakRawCount.rda"))


# merge peak
# options(scipen = 999)
merge_peak <- as.data.frame(dba_count$merged)
merge_peak$CHR <- gsub("(\\d+)","Chr\\1",merge_peak$CHR,perl = T)
merge_peak$name <- paste0(experiment_name,"_",1:dim(merge_peak)[1])
# merge_peak$scores <- "."
# merge_peak$strands <- "."
# write.table(merge_peak,
#             file = paste0("result/processed_result/", experiment_name, "_mergePeak.bed"), 
#             quote = F,sep = "\t",
#             row.names = F,col.names = F)

merge_peak_GR <- GRanges(
  seqnames = merge_peak$CHR, 
  ranges = IRanges(start = merge_peak$START, 
                   end = merge_peak$END),
  strand = "*",
  feature_id = merge_peak$name
)
assign(paste0("mergePeak_GR_", experiment_name), merge_peak_GR)

save(list = paste0("mergePeak_GR_", experiment_name),
     file = paste0("result/processed_result/", experiment_name, "_mergePeakGR.rda"))

merge_peak_GR$name <- merge_peak_GR$feature_id
rtracklayer::export.bed(con = paste0("result/processed_result/", experiment_name, "_mergePeak.bed"), 
                        object = merge_peak_GR)
