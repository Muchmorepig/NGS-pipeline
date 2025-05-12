#' @export
#' @rdnameelement
element_round_rect <- function(
    fill = NULL, colour = NULL, linewidth = NULL,
    linetype = NULL, color = NULL, inherit.blank = FALSE, radius = NULL) {
  if (!is.null(color)) colour <- color
  structure(
    list(
      fill = fill, colour = colour, linewidth = linewidth,
      linetype = linetype, radius = radius,
      inherit.blank = inherit.blank
    ),
    class = c("element_round_rect", "element_rect", "element")
  )
}

element_grob.element_round_rect <- function(
    element, x = 0.5, y = 0.5, width = 1, height = 1,
    radius = NULL,
    fill = NULL, colour = NULL, linewidth = NULL, linetype = NULL,
    ...) {
  # Thegpsettingscanoverrideelement_gp
  gp <- gpar(
    lwd = linewidth %||% element$linewidth,
    col = colour %||% element$colour,
    fill = fill %||% element$fill,
    lty = linetype %||% element$linetype
  )
  # radius
  r <- radius %||% element$radius
  roundrectGrob(x, y, width, height, r = r, gp = gp, ...)
}
