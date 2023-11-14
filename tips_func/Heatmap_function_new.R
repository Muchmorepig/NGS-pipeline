Heatmap_function <- function(matrix_data,matrix_color_fun,top_col,
                             rows_Clu = F,...){
  
  time_type <- colnames(matrix_data)
  time_type <- factor(time_type,levels = time_type)
  
  names(top_col) <- time_type
  
  
  ha_col <- HeatmapAnnotation(time = time_type,
                              col = list(
                                time = top_col
                              ),
                              gp = gpar(col = "black"),
                              show_annotation_name = F)
  
  matrix_name <- as.character(substitute(matrix_data)) 
  
  Heatmap(matrix_data,
          
          name = matrix_name,
          
          na_col = "black",
          col = matrix_color_fun,
          
          cluster_rows = rows_Clu,
          #row_split = split_type,
          border = TRUE,
          row_gap = unit(1.5, "mm"),
          row_title_rot = 0,
          row_title_gp = gpar(fontsize = 10),
          
          cluster_columns = F,
          
          ...,
          
          #show_row_names = F,
          show_column_names = F,
          
          
          top_annotation = ha_col
          #use_raster = TRUE, 
          #raster_device = "png"
          
          
          
  )
}