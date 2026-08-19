library(tidyverse)
library(patchwork)

df = read_csv("distribution_of_causal_alleles_all.csv")


df_sub = df %>%
  filter(subset_value == 0.3) %>%
  mutate(
    Architecture = factor(
      hXX,
      levels = c("h00", "h05", "h10", "h15"),
      labels = c("Recessive", "Additive", "Dominant", "Mixed")
    ),
    Population = factor(
      pop,
      levels = c("p1", "p2", "p3"),
      labels = c("African", "European", "East Asian")
    )
  )

# Summary statistics across simulations
df_summary = df_sub %>%
  group_by(Architecture, Population) %>%
  summarise(
    # Absolute counts
    mean_single   = mean(mean_count_single, na.rm = TRUE),
    sd_single     = sd(mean_count_single, na.rm = TRUE),
    se_single     = sd_single / sqrt(n()),
    mean_homo     = mean(mean_count_homo, na.rm = TRUE),
    sd_homo       = sd(mean_count_homo, na.rm = TRUE),
    se_homo       = sd_homo / sqrt(n()),
    # ROH proportions
    mean_prop_ROH = mean(prop_ROH_single, na.rm = TRUE),
    sd_prop_ROH   = sd(prop_ROH_single, na.rm = TRUE),
    se_prop_ROH   = sd_prop_ROH / sqrt(n()),
    mean_prop_ROH_homo = mean(prop_ROH_homo, na.rm = TRUE),
    sd_prop_ROH_homo   = sd(prop_ROH_homo, na.rm = TRUE),
    se_prop_ROH_homo   = sd_prop_ROH_homo / sqrt(n()),
    .groups = "drop"
  )

# FROH reference values for each population
froh_ref = tibble(
  Population = factor(
    c("African", "European", "East Asian"),
    levels = c("African", "European", "East Asian")
  ),
  FROH = c(0.061, 0.138, 0.163)
)

arch_colors = c(
  "Recessive" = '#B5134E',
  "Additive"  = '#1E88E5',
  "Dominant"  = '#FFC107',
  "Mixed"     = '#5AA99B'
)

# ---------------------------------------------------------------------------
# Panel A: Absolute number of causal alleles per individual
# ---------------------------------------------------------------------------
p_a = ggplot(df_summary, aes(x = Architecture, y = mean_single, color = Architecture)) +
  geom_pointrange(
    aes(ymin = mean_single - se_single, ymax = mean_single + se_single),
    size = 0.36, linewidth = 0.6
  ) +
  facet_wrap(~ Population) +
  scale_color_manual(values = arch_colors) +
  labs(
    x = NULL,
    y = "Mean number of causal alleles\nper individual"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position  = "none",
    strip.background = element_rect(fill = "white"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1)
  )
p_a
# ---------------------------------------------------------------------------
# Panel B: Proportion of causal alleles in ROH
# ---------------------------------------------------------------------------
p_b = ggplot(df_summary, aes(x = Architecture, y = mean_prop_ROH, color = Architecture)) +
  geom_pointrange(
    aes(ymin = mean_prop_ROH - se_prop_ROH, ymax = mean_prop_ROH + se_prop_ROH),
    size = 0.36, linewidth = 0.6
  ) +
  geom_hline(
    data = froh_ref,
    aes(yintercept = FROH),
    linetype = "dashed", color = "red", linewidth = 0.5
  ) +
  facet_wrap(~ Population) +
  scale_color_manual(values = arch_colors) +
  labs(
    x = NULL,
    y = expression(
      "Proportion of causal alleles in ROH"
    )
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position  = "none",
    strip.background = element_rect(fill = "white"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1)
  )
p_b
# ---------------------------------------------------------------------------
# Panel C: Absolute number of homozygous causal alleles per individual
# ---------------------------------------------------------------------------
p_c = ggplot(df_summary, aes(x = Architecture, y = mean_homo, color = Architecture)) +
  geom_pointrange(
    aes(ymin = mean_homo - se_homo, ymax = mean_homo + se_homo),
    size = 0.36, linewidth = 0.6
  ) +
  facet_wrap(~ Population) +
  scale_color_manual(values = arch_colors) +
  labs(
    x = NULL,
    y = "Mean number of homozygous\ncausal alleles per individual"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position  = "none",
    strip.background = element_rect(fill = "white"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1)
  )
p_c
# ---------------------------------------------------------------------------
# Panel D: Proportion of homozygous causal alleles in ROH
# ---------------------------------------------------------------------------
p_d = ggplot(df_summary, aes(x = Architecture, y = mean_prop_ROH_homo, color = Architecture)) +
  geom_pointrange(
    aes(ymin = mean_prop_ROH_homo - se_prop_ROH_homo, ymax = mean_prop_ROH_homo + se_prop_ROH_homo),
    size = 0.36, linewidth = 0.6
  ) +
  geom_hline(
    data = froh_ref,
    aes(yintercept = FROH),
    linetype = "dashed", color = "red", linewidth = 0.5
  ) +
  facet_wrap(~ Population) +
  scale_color_manual(values = arch_colors) +
  labs(
    x = NULL,
    y = "Proportion of homozygous\ncausal alleles in ROH"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position  = "none",
    strip.background = element_rect(fill = "white"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1)
  )
p_d
# ---------------------------------------------------------------------------
# Combine panels
# ---------------------------------------------------------------------------
combined = (p_a | p_b) / (p_c | p_d) +
  plot_annotation(
    title = "x=0.3",
    tag_levels = "A",
    tag_prefix = "",
    tag_suffix = "",
    theme = theme(plot.title = element_text(size = 14, hjust = 0.5, face = "bold"))
  ) &
  theme(plot.tag = element_text(face = "bold"))

combined
ggsave(
  "causal_allele_count_and_ROH_proportion.subset.0.3.pdf",
  combined,
  width = 12, height = 8
)
