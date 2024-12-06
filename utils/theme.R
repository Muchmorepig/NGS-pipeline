.mytheme <- theme_minimal() + theme(
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks = element_blank(),
  axis.title = element_blank(),
  strip.text = element_blank(),
  panel.spacing = unit(0, "lines"),
  panel.border = element_blank(),
  # panel.spacing.y = unit(-0.5, "lines"),
  panel.grid = element_blank(),
  legend.title = element_blank(),
  legend.position = "top",
  legend.background = element_rect(
    color = "#f3eaea"
  )
)
