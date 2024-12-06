make_GRange_fromMACS2peak <- function(file){
  Gr <- readPeakFile(file)
  
  seqlevels(Gr) <- gsub("Chr","",seqlevels(Gr))
  mcols(Gr) <- mcols(Gr)[1]
  colnames(mcols(Gr)) <- "feature_id"
  
  return(Gr)
}