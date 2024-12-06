plot_interest_gene <- function(interest_gene, df_select, df_select_name=""){
  
  interest_gene_cor <- df_select[df_select$gene %in% interest_gene,]
  
  if (dim(interest_gene_cor)[1] == 0){
    return(NULL)
  }  else {
    gene_cpm <- RNA.cpm[interest_gene,]
    
    peak_cpm <- data.frame(ATAC_CPM_mean[as.character(interest_gene_cor$feature_id),,
                                         drop = F],
                           check.names = F)
    peak_cpm$peak <- rownames(peak_cpm)
    peak_cpm_longer <- pivot_longer(peak_cpm,
                                    cols = -peak,
                                    names_to = "time_point",
                                    values_to = "peak_cpm")
    
    peak_cpm_RNA_cpm <- cbind(peak_cpm_longer,gene_cpm)
    peak_cpm_RNA_cpm_rvaluep <- inner_join(peak_cpm_RNA_cpm,interest_gene_cor[,c(4,6,7,8)], 
                                        by = c("peak" = "feature_id"))
    
    ggplot(data = peak_cpm_RNA_cpm_rvaluep,aes(x = gene_cpm,
                                            y = peak_cpm)) +
      geom_point(aes(fill = peak), shape = 21, size = 5, color = "black",alpha = 0.8) +
      #geom_smooth(method='lm',se = F,aes(group = peak),alpha = 0.5, color = "black") +
      stat_smooth(geom='line',method = "lm",
                  se = F,aes(group = peak),alpha = 0.5, color = "black") +
      scale_fill_manual(values = color_map) +
      theme_bw() -> p1
    
    pos_df <- group_by(peak_cpm_RNA_cpm_rvaluep,peak) %>% 
      summarise(y = mean(peak_cpm),
                x = mean(gene_cpm))
    rvalue_p_df <-  peak_cpm_RNA_cpm_rvaluep[,c(1,5,6,7)]
    
    rvalue_p_df <- rvalue_p_df[!duplicated.data.frame(rvalue_p_df), ]
    rvalue_p_df$x <- pos_df$x
    rvalue_p_df$y <- pos_df$y
    rvalue_p_df$label <- paste0("r= ",rvalue_p_df$rvalue, 
                             "\n", 
                             "p= ",rvalue_p_df$p,
                             "\n",
                             "distance= ", rvalue_p_df$distancetoTSS)
    
    p2 <- p1 + 
      geom_text_repel(aes(x=x,y=y, label=label, color=peak), data=rvalue_p_df) +
      ggtitle(label = paste0(interest_gene," ",df_select_name)) +
      scale_color_manual(values = color_map)
    
    return(p2)
  }
}
