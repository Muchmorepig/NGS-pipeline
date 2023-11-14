Heatmap_function <- function(matrix_data,matrix_color_fun,top_col,
                             rows_Clu = F,...){
  
  plant_tissue <- colnames(matrix_data)
  plant_tissue <- factor(plant_tissue,levels = plant_tissue)
  
  names(top_col) <- plant_tissue
  
  
  ha_col <- HeatmapAnnotation(plant_tissue = plant_tissue,
                              col = list(
                              plant_tissue = top_col
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
          row_title_gp = gpar(fontsize = 6),
          
          cluster_columns = F,
          
          theme(text = element_text(size = 6)),
          ...,
          
          #show_row_names = F,
          show_column_names = F,
          
          
          top_annotation = ha_col,
          use_raster = TRUE, 
          raster_device = "png"
          
          
          
  )
}