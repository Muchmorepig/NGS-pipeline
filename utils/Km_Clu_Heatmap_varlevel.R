Km_Clu_Heatmap_varlevel <- function(expr_matrix,
                           k_number,clu_save_file,
                           all_color_maps,heatmap_name,
                           level_order,
                           plot_name){
  
  # load package
  load_package <- c("BuenColors",
                    "ComplexHeatmap","circlize"
                    )
  if (!all(load_package %in% (.packages()))) lapply(load_package, library, character.only = TRUE)
  
  # Km_Clu
  ## Scale
  source("~/scripts/normal_script/normal_function/For_plot/Scale_Zscore.R")
  expr_matrix_Z <- Scale_Zscore(expr_matrix)
  
  ## Kmeans
  set.seed(19960203)
  km <- kmeans(expr_matrix_Z,centers = k_number,iter.max = 50)
  save(km,file = clu_save_file)
  
  
  # Heatmap
  
  ## set color
  names(all_color_maps) <- colnames(expr_matrix)
  
  ha_col <- HeatmapAnnotation(cell = factor(colnames(expr_matrix),
                                            levels = colnames(expr_matrix)),
                              col = list(cell = all_color_maps))
  
  ## hm
  hm <- Heatmap(expr_matrix_Z,
                       
                       column_title = heatmap_name,
                       
                       col=as.character(jdb_palette("solar_extra",type="continuous")),
                       
                       cluster_rows = FALSE, cluster_columns = FALSE, 
                       show_column_names = FALSE,show_row_names = FALSE,
                       
                       top_annotation = ha_col,
                       
                       heatmap_legend_param = list(title = "accessibilty z-socre",
                                                   title_position = "topcenter",
                                                   legend_direction = "horizontal",
                                                   border = "black",
                                                   legend_width = unit(4, "cm")),
                       
                       row_split = factor(km$cluster,levels = level_order),
                       
                       use_raster = TRUE, 
                       raster_device = "png"
  )
  
  pdf(file.path("plot",plot_name),
      width = 12, height = 12)
  draw(hm,heatmap_legend_side = "bottom")
  dev.off()
  
  
}
