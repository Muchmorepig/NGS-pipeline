# pairwise function
# also think the cofound factor

# Function to do all pair-wise DEseq comparisons
pairwise_DESeq2 <- function(x,Peak_name) {
  print(x)
  # make the contrast
  name1 <- comp[x,1] %>% as.character()
  name2 <- comp[x,2] %>% as.character()
  print(c(name1,name2,comp[x,3]))
  
  name1_condition <- sample_info$Condition[sample_info$Tissue == name1]
  name2_condition <- sample_info$Condition[sample_info$Tissue == name2]
  
  if (length(unique(name1_condition)) > 1 | 
      length(unique(name2_condition)) > 1) {
    dba_contrast <- dba.contrast(dba_count, 
                                 group1 = dba_count$masks[[name1]],
                                 group2 = dba_count$masks[[name2]],
                                 block = DBA_CONDITION,
                                 name1 = name1,
                                 name2 = name2)
    
    # diff 
    dba_diff <- dba.analyze(dba_contrast)
    
    # diff_plot
    dir.create("plot/MAplot",recursive = T)
    plot_name <- file.path("plot/MAplot",paste0(comp[x,3],".MAplot.pdf"))
    
    pdf(plot_name)
    dba.plotMA(dba_diff,fold = Fold_cutoff,
               cex.main=0.8,method=DBA_DESEQ2_BLOCK)
    abline(h = c(-Fold_cutoff,Fold_cutoff),col = "#ec008c", lty = 5)
    dev.off()
    
    # report
    dba_report_all <- dba.report(dba_diff,th = 1,
                                 method = DBA_DESEQ2_BLOCK)
    names(dba_report_all) <- paste0(Peak_name,"_",names(dba_report_all))
    dba_report_all$feature_id <- names(dba_report_all) # 加上peakID信息
  } else {
    dba_contrast <- dba.contrast(dba_count, 
                                 group1 = dba_count$masks[[name1]],
                                 group2 = dba_count$masks[[name2]],
                                 name1 = name1,
                                 name2 = name2)
    
    # diff 
    dba_diff <- dba.analyze(dba_contrast)
    
    # diff_plot
    dir.create("plot/MAplot",recursive = T)
    plot_name <- file.path("plot/MAplot",paste0(comp[x,3],".MAplot.pdf"))
    
    pdf(plot_name)
    dba.plotMA(dba_diff,fold = Fold_cutoff,
               cex.main=0.8)
    abline(h = c(-Fold_cutoff,Fold_cutoff),col = "#ec008c", lty = 5)
    dev.off()
    
    # report
    dba_report_all <- dba.report(dba_diff,th = 1)
    names(dba_report_all) <- paste0(Peak_name,"_",names(dba_report_all))
    dba_report_all$feature_id <- names(dba_report_all) # 加上peakID信息
    
  }

  # 把seqlevel改成Ensemble的格式
  seqlevels(dba_report_all) <- gsub("Chr","",seqlevels(dba_report_all))
  
  # annotate
  peakAnno <- annotatePeak(dba_report_all,
                           TxDb = TxDb.Athaliana.BioMart.plantsmart28,
                           level = "gene")
  # add geneSymbol
  peakAnno_geneSymbol <- left_join(as.data.frame(peakAnno),gene_alias)
  peakAnno_geneSymbol$seqnames <- gsub("(\\d+)","Chr\\1",peakAnno_geneSymbol$seqnames)
  
  # result output
  dir.create("result/Diff_Peak_Anno",recursive = T)
  result_name <- file.path("result/Diff_Peak_Anno",
                           paste0(experiment_name,"_",comp[x,3],".csv")
  )
  
  # 最好quote是TRUE，因为有些是分隔符分隔的
  write.csv(peakAnno_geneSymbol,result_name,row.names = F)
  
  # compare GO analysis
  peakAnno_diff_gene_list <- list()
  
  peakAnno_diff_gene_list[["up"]] <- subset(peakAnno@anno,Fold > Fold_cutoff & FDR < FDR_cutoff)$geneId %>% unique()
  peakAnno_diff_gene_list[["down"]] <- subset(peakAnno@anno,Fold < -Fold_cutoff & FDR < FDR_cutoff)$geneId %>% unique()
  
  
  cc <- try(compareCluster(geneClusters = peakAnno_diff_gene_list,
                           fun="enrichGO", OrgDb = org.At.tair.db,
                           keyType = "TAIR",
                           ont = "BP",pvalueCutoff=0.05), # 阈值卡的大点
            silent = T)
  if ('try-error' %in% class(cc)) {
    print("No_cc")
  } else {
    ## ouputname
    ## ouputname
    dir.create("result/Diff_Peak_GO_result",recursive = T)
    dir.create("plot/Diff_Peak_GO_plot",recursive = T)
    compareGO_plot_outputname <- file.path("plot/Diff_Peak_GO_plot",paste0(experiment_name,"_",comp[x,3],"_compareGO.pdf"))
    compareGO_result_outputname <- file.path("result/Diff_Peak_GO_result",paste0(experiment_name,"_",comp[x,3],"_compareGO.csv"))
    
    
    ## output
    write.csv(cc@compareClusterResult %>% as.data.frame(),
              compareGO_result_outputname,row.names = F)
    
    p <- dotplot(cc, showCategory=20)
    
    pdf(compareGO_plot_outputname,width = 10,height = 10)
    print(p)
    dev.off()
    
  }
  
}