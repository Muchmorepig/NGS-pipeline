library(grid)
word_cloud_grob = function(text, fontsize,
                           line_space = unit(4, "pt"), word_space = unit(4, "pt"), max_width = unit(80, "mm"),
                           col = function(fs) circlize::rand_color(length(fs), luminosity = "dark"),
                           test = FALSE) { # width in mm
  
  if(length(text) != length(fontsize)) {
    stop("`text` and `fontsize` should the same length.")
  }
  
  od = order(fontsize, decreasing = TRUE)
  text = text[od]
  fontsize = fontsize[od]
  
  if(Sys.info()["sysname"] == "Darwin" && dev.interactive()) {
    ComplexHeatmap:::dev.null()
    on.exit(ComplexHeatmap:::dev.off2())
  }
  
  n = length(text)
  text_gb_lt = lapply(seq_len(n), function(i) textGrob(text[i], gp = gpar(fontsize = fontsize[i])))
  text_width = vapply(text_gb_lt, function(gb) convertWidth(grobWidth(gb), "mm", valueOnly = TRUE), 0)
  text_height = vapply(text_gb_lt, function(gb) convertHeight(grobHeight(gb), "mm", valueOnly = TRUE), 0)
  
  if(is.unit(line_space)) line_space = convertHeight(line_space, "mm", valueOnly = TRUE)
  if(is.unit(word_space)) word_space = convertWidth(word_space, "mm", valueOnly = TRUE)
  
  x = numeric(n)
  y = numeric(n)
  current_line_height = 0
  current_line_width = 0
  
  # the first text
  current_line_height = text_height[1]
  current_line_width = text_width[1]
  x[1] = 0
  y[1] = 0
  
  w = text_width[1]
  h = text_height[1]
  
  if(is.unit(max_width)) {
    max_width = convertWidth(max_width, "mm", valueOnly = TRUE)
  }
  
  for(i in seq_len(n)[-1]) {
    # the next text can be put on the same line
    if(current_line_width + text_width[i] <= max_width) {
      x[i] = current_line_width + word_space
      y[i] = y[i-1] # same as previous one
      current_line_width = x[i] + text_width[i]
      w = max(w, current_line_width)
      h = max(h, y[i] + text_height[i])
    } else { # the next text need to be put on the next line
      x[i] = 0
      y[i] = current_line_height + line_space
      current_line_width = text_width[i]
      current_line_height = y[i] + text_height[i]
      w = max(w, current_line_width)
      h = max(h, current_line_height)
    }
  }
  
  if(is.character(col) || is.numeric(col)) {
    if(length(col) == 1) col = rep(col, n)
    col_fun = function(fontsize) return(col)
  } else if(is.function(col)) {
    col_fun = col
  } else {
    stop("`col` can only be a function or a character vector.")
  }
  
  if(test) {
    gl = gList(
      rectGrob(),
      textGrob(text, x = x, y = y, gp = gpar(fontsize = fontsize, col = col_fun(fontsize)),
               default.units = "mm", just = c(0, 0)),
      rectGrob(x = x, y = y, width = text_width, height = text_height, default.units = "mm", just = c(0, 0))
      
    )
  } else {
    gl = gList(
      textGrob(text, x = x, y = y, gp = gpar(fontsize = fontsize, col = col_fun(fontsize)),
               default.units = "mm", just = c(0, 0))
    )
  }
  
  gb = gTree(children = gl, cl = "word_cloud", vp = viewport(width = unit(w, "mm"), height = unit(h, "mm")))
  return(gb)
}

widthDetails.word_cloud = function(x) {
  x$vp$width
}

heightDetails.word_cloud = function(x) {
  x$vp$height
}

scale_fontsize = function(x, rg = c(1, 30), fs = c(4, 16)) {
  k = (fs[2] - fs[1])/(rg[2] - rg[1])
  b = fs[2] - k*rg[2]
  y = k*x + b
  y[y < fs[1]] = fs[1]
  y[y > fs[2]] = fs[2]
  round(y)
}