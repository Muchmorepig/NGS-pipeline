one_one <- function(x,
                    Type = "type",
                    cutoff_padj = 0.05,
                    comp = comp_info,
                    dds_choose = dds,
                    od = outdir,
                    shrink = TRUE,
                    seq = "ATAC") {
  # Do the DESeq contrast
  require(DESeq2)
  require(dplyr)

  treat <- x[1]
  control <- x[2]
  key <- x[3]
  message(c(treat, "\t", control, "\t", key))

  suppressMessages(
    if (shrink) {
      res_lfc <- lfcShrink(
        dds = dds_choose,
        type = "ashr",
        alpha = cutoff_padj,
        contrast = c(Type, treat, control)
      )
    } else {
      res_lfc <- results(
        object = dds_choose,
        alpha = cutoff_padj,
        contrast = c(Type, treat, control)
      )
    }
  )

  if (seq == "ATAC") {
    res_lfc %>%
      as_tibble(rownames = "feature_id") %>%
      mutate(VS = key) %>%
      dplyr::select(1, 3, 5:7) %>%
      readr::write_csv(file = file.path(od, paste0(key, "_DA.csv")))
  } else if (seq == "RNA") {
    res_lfc %>%
      as_tibble(rownames = "geneId") %>%
      mutate(VS = key) %>%
      dplyr::select(1, 3, 5:7) %>%
      readr::write_csv(file = file.path(od, paste0(key, "_DE.csv")))
  }

  return(" ^_^ ")
}
