PCAplot <- function(object, intgroup = c("group"), 
                    ntop = 2000, size = 5, 
                    pc = c("PC1","PC2"), 
                    returnData = F,
                    show = T,
                    LP = "top"){
  
  if(!"ggrepel" %in% (.packages())) library(ggrepel)
  
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = T)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  names(percentVar) <- colnames(pca$x)
  
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[,intgroup, drop = F])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = ":"))
    }else{
    colData(object)[[intgroup]]
      }

  if (length(pc) == 2){
    
    d <- data.frame(PCa = pca$x[, pc[1]], PCb = pca$x[, pc[2]], 
                    group = group, 
                    intgroup.df, 
                    name = colnames(object))
    
    p <- ggplot(data = d, aes_string(x = "PCa", y = "PCb", color = "group")) + 
      geom_point(size = size) + 
      xlab(paste0(pc[1], ": ", round(percentVar[pc[1]] * 100), "% variance")) + 
      ylab(paste0(pc[2], ": ", round(percentVar[pc[2]] * 100), "% variance")) + 
      coord_fixed() +
      theme_bw() +
      theme(legend.background = element_blank(),
            legend.position = LP) +
      geom_label_repel(aes(label=name), 
                       fontface="bold", 
                       color="grey50", 
                       box.padding=unit(0.35, "lines"), 
                       point.padding=unit(0.5, "lines"),
                       segment.colour = "grey50",
                       max.overlaps = 15)+
      guides(color = guide_legend(title = "Group"))
    if (show){
      print(p)
    }
    if (returnData) {
      #attr(d, "percentVar") <- percentVar[pc[1:length(pc)]]
      return(
        list(PCA  = pca,
             pData= d,
             plot = p,
             percentVar = percentVar)
      )}
    
  }else{
    if(!"plotly" %in% (.packages())) library(plotly)
    d <- data.frame(PCa = pca$x[, pc[1]], PCb = pca$x[, pc[2]], PCc = pca$x[, pc[3]],
                    group = group, 
                    intgroup.df, 
                    name = colnames(object))
    colnames(d)[1:length(pc)] <- pc
    
    d %>% 
      plot_ly(x = ~PC1, y = ~PC2, z = ~PC3) %>%
      add_markers(
        color = ~group,) %>% 
      add_text(x = ~PC1, y = ~PC2, z = ~PC3, text = ~name) %>% 
      layout(
        scene = list(
          xaxis = list(title = paste0("PC1: ",round(percentVar["PC1"] * 100), "% variance")),
          yaxis = list(title = paste0("PC2: ",round(percentVar["PC2"] * 100), "% variance")),
          zaxis = list(title = paste0("PC3: ",round(percentVar["PC3"] * 100), "% variance"))
        )
      ) -> p 
    if (show){
      print(p)
    }
    if (returnData) {
      #attr(d, "percentVar") <- percentVar[pc[1:length(pc)]] 
      return(list(PCA  = pca,
                  pData= d,
                  plot = p,
                  percentVar = percentVar))
    }
  }
}
