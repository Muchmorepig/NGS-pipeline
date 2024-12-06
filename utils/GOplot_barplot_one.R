GOplot_barplot <- function(data, mapping = aes(
                                 x = reorder(Description, -log10(p.adjust)),
                                 y = -log10(p.adjust),
                                 label = Description
                           ),
                           barfill = "#ff4d79", baralpha = 0.9,
                           title_name = "") {
      if (!"BuenColors" %in% (.packages())) library(BuenColors)

      ggplot(data, mapping) +
            geom_bar(
                  stat = "identity",
                  position = position_stack(reverse = T),
                  fill = barfill, alpha = baralpha
            ) +
            theme_bw() +
            pretty_plot(fontsize = 12) +
            coord_flip() +
            L_border() +
            labs(x = "", y = "-log10(padj)", title = title_name) +
            theme(
                  legend.position = "none",
                  axis.text.y = element_blank(),
                  axis.ticks.y = element_blank()
            ) +
            geom_text(
                  color = "black",
                  position = position_stack(vjust = 0), hjust = 0, size = 5
            ) +
            scale_y_continuous(expand = c(0, 0.1))
}
