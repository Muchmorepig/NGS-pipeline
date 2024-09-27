#' @export
one_one <- function(x,
                    Type = "type",
                    cutoff_padj = 0.05,
                    comp = comp_info,
                    dds_choose = dds,
                    od = outdir,
                    shrink = TRUE,
                    seq = "ATAC") {
  require(DESeq2)
  require(dplyr)

  treat <- x[1]
  control <- x[2]
  key <- x[3]
  message(c(treat, "\t", control, "\t", key))

  res_lfc <- suppressMessages({
    if (shrink) {
      lfcShrink(
        dds = dds_choose,
        type = "ashr",
        alpha = cutoff_padj,
        contrast = c(Type, treat, control)
      )
    } else {
      results(
        object = dds_choose,
        alpha = cutoff_padj,
        contrast = c(Type, treat, control)
      )
    }
  })

  if (is.null(res_lfc)) {
    warning("Results are NULL for comparison: ", key)
    return(NULL)
  }

  output_filename <- switch(seq,
    "ATAC" = paste0(key, "_DA.csv"),
    "RNA" = paste0(key, "_DE.csv"),
    paste0(key, "_results.csv")
  ) # Default or fallback file name

  res_lfc %>%
    as_tibble(rownames = ifelse(seq == "ATAC", "feature_id", "geneId")) %>%
    mutate(VS = key) %>%
    dplyr::select(1, 3, 5:7) %>%
    readr::write_csv(file = file.path(od, output_filename))

  message("File written: ", output_filename)
}
