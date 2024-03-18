file.move <- function(files, new_path_prefix) {
  for (i in seq_along(files)) {
    file_name <- basename(files[1])
    new_path <- file.path(new_path_prefix, file_name)

    file.link(files[i], new_path)
    unlink(files[i])

    files[i] <- new_path
  }

  return(files)
}
