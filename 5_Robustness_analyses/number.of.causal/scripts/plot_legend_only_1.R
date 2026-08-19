library(tidyverse)
library(cowplot)

arch_colors <- c(
  "Recessive" = '#B5134E',
  "Additive"  = '#1E88E5',
  "Dominant"  = '#FFC107',
  "Mixed"     = '#5AA99B'
)

dummy <- tibble(
  x = 1:4,
  y = 1:4,
  Architecture = factor(
    c("Recessive", "Mixed", "Additive", "Dominant"),
    levels = c("Recessive", "Mixed", "Additive", "Dominant")
  ),
  FROH = 0.5
)

p <- ggplot(dummy) +
  geom_pointrange(
    aes(x = x, y = y, ymin = y - 0.1, ymax = y + 0.1, color = Architecture),
    size = 0.5, linewidth = 0.6
  ) +
  geom_hline(aes(yintercept = FROH, linetype = "FROH"), color = "red", linewidth = 0.5) +
  scale_color_manual(values = arch_colors) +
  scale_linetype_manual(
    name   = NULL,
    values = c("FROH" = "dashed"),
    labels = expression(F[ROH])
  ) +
  guides(
    color    = guide_legend(order = 1),
    linetype = guide_legend(order = 2, override.aes = list(color = "red"))
  ) +
  theme(
    legend.text       = element_text(size = 11),
    legend.title      = element_text(size = 12),
    legend.key.width  = unit(1.5, "cm"),
    legend.key        = element_rect(fill = "transparent", colour = NA),
    legend.background = element_rect(fill = "transparent", colour = NA)
  )

legend <- get_legend(p)

ggsave(
  "legend_only.pdf",
  plot_grid(legend),
  width = 6, height = 2,
  bg = "transparent"
)
