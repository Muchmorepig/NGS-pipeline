#' Select Color Palette
#'
#' This function selects a color palette based on the specified palette name.
#' @param palette A character string indicating the color palette to use.
#' Possible values are "col29", "col28", "col_sim", "col_sex", and "col_33".
#' @return A vector of colors corresponding to the selected palette.
#' @export
select_colors <- function(palette) {
  switch(palette,
    col29 = c(
      "#eb9aeb", "#0d9eff", "#e8b964", "#465ca8", "#089ac3", "#dacd50",
      "#0b6a24", "#66c4cc", "#4bace0", "#907bb6", "#5bb449", "#c5da6b",
      "#ef9d3b", "#86c69e", "#6f6f70", "#68ac7d", "#557d65", "#7c4ca2",
      "#f1d1a0", "#649da8", "#ff9ea9", "#c3e1a2", "#c779b7", "#727e57",
      "#a791a8", "#ae9787", "#008c99", "#f46952", "#c94e91"
    ),
    col28 = c(
      "#d60101", "#d0595a", "#f786a8", "#f74a3a",
      "#faad89", "#f77c43", "#f85408", "#f01778",
      "#b7889d", "#ecd0d0", "#17e7ee", "#8dcaec",
      "#04a5fc", "#66c2a5", "#3bd608", "#a9e097",
      "#2c6917", "#cab9f5", "#df8633", "#05968e",
      "#5e6166", "#1872a3", "#7774c2", "#393b78",
      "#0b13f1", "#a00b98", "#63065d", "#2e012e"
    ),
    col_6 = c(
      "#49a0b3", "#a0d9a4", "#fee998", "#ffb763", "#f1949c", "#cd281c"
    ),
    col_2 = c("#fd81ae", "#4fd5fa"),
    col_33 = c(
      "#9ec4d3", "#16c1c7", "#f786a8", "#faad89", "#fdd0a2",
      "#f77c43", "#f85408", "#05968e", "#b7889d", "#ecd0d0",
      "#c6dbef", "#a9e097", "#04a5fc", "#66c2a5", "#74c476",
      "#8dcaec", "#f01778", "#cab9f5", "#df8633", "#2c6917",
      "#d0595a", "#1872a3", "#7774c2", "#a00b98", "#fa6ac3",
      "#d29c9d", "#cf3c29", "#d2d2e4", "#63065d", "#2e012e",
      "#701a12", "#393b78", "#363cea"
    ),
    stop("Palette not recognized. Available options: col29, col28, col_sim, col_sex, col_33")
  )
}
