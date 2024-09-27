#' @export
emb_top_middle_low_cells <- function(obj, layer_name, sel_genes, ps = 0.3, return_data = FALSE) {
  v <- stringr::str_sub(as.character(obj@version), 1, 1)

  message("Seurat Version: ", obj@version)

  if (v == "4") {
    mtx <- GetAssayData(obj, layer_name)
  } else {
    mtx <- LayerData(obj, layer_name)
  }

  sub_mtx <- mtx[rownames(mtx) %in% sel_genes, ]

  temp_mtx <- colSums(sub_mtx) / colSums(mtx)
  temp_mtx <- temp_mtx[temp_mtx > 0]

  len <- length(temp_mtx) / 5

  sorted_temp_mtx <- sort(-temp_mtx)
  top_cells <- names(sorted_temp_mtx[1:round(len / 3)])
  middle_cells <- names(sorted_temp_mtx[(round(len / 3) + 1):(round(len / 3) * 2)])
  low_cells <- names(sorted_temp_mtx[((round(len / 3) * 2) + 1):len])

  obj[["ident_1"]] <- ""

  obj$ident_1[colnames(obj) %in% low_cells] <- "1_low"
  obj$ident_1[colnames(obj) %in% middle_cells] <- "2_middle"
  obj$ident_1[colnames(obj) %in% top_cells] <- "3_top"

  p <- DimPlot(
    object = obj,
    reduction = "umap",
    group.by = "ident_1",
    repel = TRUE,
    raster = FALSE,
    pt.size = ps,
    cols = c("lightgrey", "orangered1", "orangered3", "orangered4"),
    order = TRUE
  ) + theme(plot.title = element_blank())
  
  if (return_data) {
    return(list(plot = p, data = obj$ident_1))
  } else {
    return(p)
  }
}
