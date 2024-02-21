library(ArchR)

proj_path <- "./Result/Merge/"

proj <- loadArchRProject(proj_path)

colnames(proj@cellColData)

getGroupBW(
  ArchRProj = proj,
  groupBy = "ClusterHarmony_iter4",
  normMethod = "ReadsInTSS",
)
