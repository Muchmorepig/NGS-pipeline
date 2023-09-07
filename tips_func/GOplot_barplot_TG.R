GOplot_barplot_group <- function(data,
                                 color_map = c("#256b99","#db2029"),
                                 baralpha = 0.85,
                                 title_name = ""){
  
  # if (!"BuenColors" %in% (.packages())) library(BuenColors)
  if (!"dplyr" %in% (.packages())) library(dplyr)
  
  data <- data %>% arrange(.,Cluster,-log10(p.adjust))
  data$Description <- factor(data$Description, levels = data$Description)
  
  index <- table(data$Cluster) %>% as.numeric() %>% "["(1)
  
  data$x_pos <- as.numeric(data$Description)
  data$x_pos <- ifelse(data$x_pos > index, data$x_pos + 0.5, data$x_pos)
  
  
  ggplot(data,aes(x = x_pos, y = -log10(p.adjust),
                  label = Description,
                  fill = Cluster)) +
    geom_bar(stat = "identity",
             position = position_dodge(width = 8),
             alpha = 0.85)+
    coord_flip() +
    L_border() +
    labs(x="",y="-log10(padj)",title = title_name) +
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
    guides(fill = guide_legend(reverse = T)) +
    geom_text(color="black",
              position = position_stack(vjust = 0),hjust=0,size=3.5) +
    scale_y_continuous(expand = c(0,0.1)) +
    scale_x_continuous(expand = c(0,0.2)) +
    scale_fill_manual(values = color_map)
  
  

    
}

  
  
          
