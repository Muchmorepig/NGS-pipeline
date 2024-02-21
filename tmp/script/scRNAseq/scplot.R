features = c("AT3G24140","AT3G26744","AT2G46720","AT1G80080",
             "AT4G38160","AT1G10070","AT1G52270","AT4G18290",
             "AT4G24060","AT5G54510","AT2G30750","AT4G21750",
             "ATCG00340","AT1G67090","AT3G41768","AT1G29910",
             "AT1G67090","AT1G29910","AT1G29920","AT1G29930",
             "AT1G15820","AT5G65590")
outdir <- "pngs"
dir.create(outdir)
for ( feature in features){
  
  fn <- file.path(outdir, paste0(feature, ".png"))
  
  png(fn, width = 1000, height = 800)
  
  p1 <- VlnPlot(seurat.obj,pt.size = 0,split.by = "orig.ident", features = feature)
  p2 <- FeaturePlot(seurat.obj, split.by = "orig.ident", features = feature)
  print(p1 + p2)
  dev.off()
}
