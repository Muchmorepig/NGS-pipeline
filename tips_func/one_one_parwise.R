one_one_parwise <- function(x, Type, dds_choose) {
  print(x)
  # Do the DESeq contrast
  treat <- as.character(compare_expand[x,1])
  control <- as.character(compare_expand[x,2])
  key <- as.character(compare_expand[x,3])
  print(c(treat, control, key))
  
  res <- results(dds_choose, 
                 contrast = c(Type, treat, control), 
                 alpha = cutoff_padj)
  
  
  res_lfc <- lfcShrink(dds = dds_choose,
                       res = res,
                       type = "ashr")
  
  # output res result
  count_norm_res <- count_norm[, c(which(data_type == treat),
                                   which(data_type == control))] %>% 
    as_tibble(rownames = "feature_id")
  
  count_norm_res %>% 
    inner_join(as_tibble(res_lfc, rownames = "feature_id")) %>% 
    dplyr::select(-c("baseMean", "lfcSE")) -> DESeq2_result
  
  # MA-plot
  DESeq2_result %>% 
    dplyr::select(2:6, 8) %>% 
    mutate(mean = rowMeans(dplyr::select(., 1:4))) %>% 
    mutate(sig = case_when(
      abs(log2FoldChange) > cutoff_lfc & padj < cutoff_padj ~ TRUE,
      TRUE ~ FALSE
    )) %>% 
    dplyr::select(7, 5, 8) -> plotMA_data
  
  plotMA_data %>% 
    count(sig) %>% 
    pull(n) -> summarise_data
  
  if (length(summarise_data) == 1) {
    summarise_data[2] <- 0
  }
  
  percent <- (round(summarise_data[2] / sum(summarise_data), digits = 3)) * 100
  
  pdf(paste0("plot/MAplot/", key, ".pdf"))
  geneplotter::plotMA(plotMA_data)
  title(paste0("total peak:", sum(summarise_data), "    ",
               "sig peak:", summarise_data[2] , "(", percent, "%)", "\n",
               "Fold cutoff: ", cutoff_lfc, "    pdj cutoff: ", cutoff_padj))
  dev.off()
    
  
  peakAnno %>% 
    inner_join(DESeq2_result) %>% 
    dplyr::select(1:4, 21:27, 6:20) -> diff_peak_anno
  
  readr::write_csv(diff_peak_anno, file = paste0("result/Diff_Peak_Anno/", 
                                                 key, 
                                                 "_diff_result.csv"))
  
  
  # GO result
  diff_gene <- list()
  diff_gene[["up"]] <- diff_peak_anno %>% 
    filter(padj < cutoff_padj, 
           log2FoldChange > cutoff_lfc) %>% 
    pull(geneId) %>% 
    unique()
    
  diff_gene[["down"]] <- diff_peak_anno %>% 
    filter(padj < cutoff_padj, 
           log2FoldChange < -cutoff_lfc) %>% 
    pull(geneId) %>% 
    unique()
  
  
  cc <- try(clusterProfiler::compareCluster(geneClusters = diff_gene,
                           fun = "enrichGO",
                           OrgDb = org.At.tair.db,
                           keyType = "TAIR",
                           ont = "BP"),
            silent = T)
  if ('try-error' %in% class(cc)) {
    print("No_cc")
  } else {
    readr::write_csv(cc@compareClusterResult, 
                     file = paste0("result/Diff_Peak_GO_result/", key, 
                                "_GO_result.csv"))
  }
}
