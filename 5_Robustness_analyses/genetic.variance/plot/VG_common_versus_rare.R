library(tidyverse)

# x = 0.3 -----------------------------------------------------------------
x = '0.3'
df = read.csv(paste0("all_simulations_variance_partition_common_versus_rare.",x,".csv"), check.names = FALSE)
df_mean <- df %>%
  select(-`Sum_Check (%)`) %>%
  mutate(
    prop_less_than_0_001_V_G = less_than_0_001_V_G / Total_V_G,
    prop_between_0_001_and_0_01_V_G = between_0_001_and_0_01_V_G / Total_V_G,
    prop_more_than_0_01_V_G = more_than_0_01_V_G / Total_V_G,
    prop_Covariance = Covariance / Total_V_G
  ) %>%
  select(-less_than_0_001_V_G, -between_0_001_and_0_01_V_G, -more_than_0_01_V_G, -Covariance, -Total_V_G) %>%
  group_by(hXX, pop, tau) %>%
  summarise(
    mean_prop_less_than_0_001_V_G  = 100*mean(prop_less_than_0_001_V_G, na.rm = TRUE),
    mean_prop_between_0_001_and_0_01_V_G = 100*mean(prop_between_0_001_and_0_01_V_G, na.rm = TRUE),
    mean_prop_more_than_0_01_V_G = 100*mean(prop_more_than_0_01_V_G, na.rm = TRUE),
    mean_prop_Covariance = 100*mean(prop_Covariance, na.rm = TRUE),
    .groups       = "drop"
  )%>%
  mutate(
    hXX = case_when(
      hXX == "h00" ~ "Recessive",
      hXX == "h05" ~ "Additive",
      hXX == "h10" ~ "Dominant",
      hXX == 'h15' ~ "Mixed",
    ),
    hXX = factor(hXX, levels = c("Recessive", "Additive", "Dominant", "Mixed"))
  )

df_long <- df_mean %>%
  pivot_longer(
    cols = c(mean_prop_less_than_0_001_V_G, mean_prop_between_0_001_and_0_01_V_G, mean_prop_more_than_0_01_V_G, mean_prop_Covariance),
    names_to = "Component",
    values_to = "Percentage"
  ) %>%
  mutate(
    Component = factor(Component, levels = c("mean_prop_less_than_0_001_V_G", "mean_prop_between_0_001_and_0_01_V_G", "mean_prop_more_than_0_01_V_G", "mean_prop_Covariance"))
  )

p <- ggplot(df_long, aes(x = factor(tau), y = Percentage, fill = Component)) +
  geom_col(position = position_stack(reverse = TRUE), width = 0.5, color = "black", linewidth = 0.2) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.6) +
  geom_hline(yintercept = 100, linetype = "dashed", color = "black", linewidth = 0.6, alpha = 0.7) +
  
  facet_grid(hXX ~ factor(pop, levels = c("African", "European", "East Asian"))) +
  
  scale_fill_manual(
    values = c(
      "mean_prop_less_than_0_001_V_G" = "#56B4E9", 
      "mean_prop_between_0_001_and_0_01_V_G"   = "#CC79A7",  
      "mean_prop_more_than_0_01_V_G" = "#E69F00",
      "mean_prop_Covariance" = "#808080"
    ),
    labels = c(
      "mean_prop_less_than_0_001_V_G" = "f < 0.1%",
      "mean_prop_between_0_001_and_0_01_V_G"   = "0.1% \u2264 f < 1%",
      "mean_prop_more_than_0_01_V_G" = "1% \u2264 f",
      "mean_prop_Covariance" = "Covariance (LD)"
    )
  ) +
  
  labs(
    x = "Tau",
    y = expression(bold("Proportion of V"[G]~"explained")),
    title = paste0("x = ",x),
    fill = "Variance Component"
  ) +
  
  theme_bw(base_size = 12.5, base_family = "Arial") +
  theme(
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text = element_text(face = "bold", size = 12),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(face = "bold", size = 15),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16)
  )
p
ggsave(paste0("Genetic_Variance_Rare_versus_Common.",x,".pdf"), plot = p, width = 8, height = 8, device = cairo_pdf)


# x = 1.0 -----------------------------------------------------------------
x = '1.0'
df = read.csv(paste0("all_simulations_variance_partition_common_versus_rare.",x,".csv"), check.names = FALSE)
df_mean <- df %>%
  select(-`Sum_Check (%)`) %>%
  mutate(
    prop_less_than_0_001_V_G = less_than_0_001_V_G / Total_V_G,
    prop_between_0_001_and_0_01_V_G = between_0_001_and_0_01_V_G / Total_V_G,
    prop_more_than_0_01_V_G = more_than_0_01_V_G / Total_V_G,
    prop_Covariance = Covariance / Total_V_G
  ) %>%
  select(-less_than_0_001_V_G, -between_0_001_and_0_01_V_G, -more_than_0_01_V_G, -Covariance, -Total_V_G) %>%
  group_by(hXX, pop, tau) %>%
  summarise(
    mean_prop_less_than_0_001_V_G  = 100*mean(prop_less_than_0_001_V_G, na.rm = TRUE),
    mean_prop_between_0_001_and_0_01_V_G = 100*mean(prop_between_0_001_and_0_01_V_G, na.rm = TRUE),
    mean_prop_more_than_0_01_V_G = 100*mean(prop_more_than_0_01_V_G, na.rm = TRUE),
    mean_prop_Covariance = 100*mean(prop_Covariance, na.rm = TRUE),
    .groups       = "drop"
  )%>%
  mutate(
    hXX = case_when(
      hXX == "h00" ~ "Recessive",
      hXX == "h05" ~ "Additive",
      hXX == "h10" ~ "Dominant",
      hXX == 'h15' ~ "Mixed",
    ),
    hXX = factor(hXX, levels = c("Recessive", "Additive", "Dominant", "Mixed"))
  )

df_long <- df_mean %>%
  pivot_longer(
    cols = c(mean_prop_less_than_0_001_V_G, mean_prop_between_0_001_and_0_01_V_G, mean_prop_more_than_0_01_V_G, mean_prop_Covariance),
    names_to = "Component",
    values_to = "Percentage"
  ) %>%
  mutate(
    Component = factor(Component, levels = c("mean_prop_less_than_0_001_V_G", "mean_prop_between_0_001_and_0_01_V_G", "mean_prop_more_than_0_01_V_G", "mean_prop_Covariance"))
  )

p <- ggplot(df_long, aes(x = factor(tau), y = Percentage, fill = Component)) +
  geom_col(position = position_stack(reverse = TRUE), width = 0.5, color = "black", linewidth = 0.2) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.6) +
  geom_hline(yintercept = 100, linetype = "dashed", color = "black", linewidth = 0.6, alpha = 0.7) +
  
  facet_grid(hXX ~ factor(pop, levels = c("African", "European", "East Asian"))) +
  
  scale_fill_manual(
    values = c(
      "mean_prop_less_than_0_001_V_G" = "#56B4E9", 
      "mean_prop_between_0_001_and_0_01_V_G"   = "#CC79A7",  
      "mean_prop_more_than_0_01_V_G" = "#E69F00",
      "mean_prop_Covariance" = "#808080"
    ),
    labels = c(
      "mean_prop_less_than_0_001_V_G" = "f < 0.1%",
      "mean_prop_between_0_001_and_0_01_V_G"   = "0.1% \u2264 f < 1%",
      "mean_prop_more_than_0_01_V_G" = "1% \u2264 f",
      "mean_prop_Covariance" = "Covariance (LD)"
    )
  ) +
  
  labs(
    x = "Tau",
    y = expression(bold("Proportion of V"[G]~"explained")),
    title = paste0("x = ",x),
    fill = "Variance Component"
  ) +
  
  theme_bw(base_size = 12.5, base_family = "Arial") +
  theme(
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text = element_text(face = "bold", size = 12),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(face = "bold", size = 15),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16)
  )
p
ggsave(paste0("Genetic_Variance_Rare_versus_Common.",x,".pdf"), plot = p, width = 8, height = 8, device = cairo_pdf)
