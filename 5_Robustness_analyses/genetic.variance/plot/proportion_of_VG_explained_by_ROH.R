library(dplyr)
library(tidyr)
library(ggplot2)

# x=0.3(centered and normalized) ------------------------------------------
df_centered_normalized = read.csv("ALL_variance_partition_summary_centered_normalization.csv")
df_centered_normalized_0.3 = df_centered_normalized %>%
  filter(subset_value == 0.3) %>% 
  mutate(
    Prop_VA = Total_V_A / Total_V_G,
    Prop_VD = Total_V_D / Total_V_G,
    Prop_Cov = Total_cov / Total_V_G,
    
    Prop_ROH_VG = ROH_V_G / Total_V_G,
    Prop_NonROH_VG = Non_ROH_V_G / Total_V_G,
    Prop_ROH_Cov = 2*ROH_cov / Total_V_G
  ) %>%
  group_by(hXX, pop, tau, subset_value) %>%
  summarise(across(-simID, ~mean(.x, na.rm = TRUE)), .groups = 'drop') %>%
  select(-Total_check, -ROH_check) %>%
  mutate(
    hXX = case_when(
      hXX == "h00" ~ "Recessive",
      hXX == "h05" ~ "Additive",
      hXX == "h10" ~ "Dominant",
      hXX == 'h15' ~ "Mixed",
    ),
    hXX = factor(hXX, levels = c("Recessive", "Additive", "Dominant", "Mixed"))
  )%>%
  mutate(
    pop = case_when(
      pop == "p1" ~ "African",
      pop == "p2" ~ "European",
      pop == "p3" ~ "East Asian",
    ),
    pop = factor(pop, levels = c("African", "European", "East Asian"))
  )

df_centered_normalized_0.3_long = df_centered_normalized_0.3 %>%
  select(hXX, pop, tau, Prop_ROH_VG, Prop_NonROH_VG, Prop_ROH_Cov) %>%
  pivot_longer(
    cols = c(Prop_ROH_VG, Prop_NonROH_VG, Prop_ROH_Cov),
    names_to = "region",
    values_to = "proportion"
  )%>%
  mutate(region = factor(region, levels = c("Prop_ROH_VG", "Prop_NonROH_VG", "Prop_ROH_Cov")))

plot_centered_normalized_0.3 = ggplot(df_centered_normalized_0.3_long, aes(x = factor(tau), y = proportion, fill = region)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(hXX ~ pop) +
  scale_fill_manual(values = c(
    "Prop_ROH_VG" = "#A9DEF9",      
    "Prop_NonROH_VG" = "#90A955",   
    "Prop_ROH_Cov" = "#5E6472"     
  ),
  labels = c(
    "Prop_ROH_VG" = "ROH",
    "Prop_NonROH_VG" = "Non-ROH",
    "Prop_ROH_Cov" = "Covariance"
  )) +
  labs(
    x = "Tau",
    y = expression(bold("Proportion of V"[G]~"explained")),
    fill = "Region",
    title = "x = 0.3"
  ) +
  coord_cartesian(ylim = c(-0.05, 1.05)) +
  theme_bw() +
  theme(
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)) 
plot_centered_normalized_0.3
ggsave("Genetic_Variance_explained_by_ROH_centered_normalized_0.3.pdf", plot = plot_centered_normalized_0.3, width = 8, height = 8, device = cairo_pdf)
# x=1.0(centered and normalized) ------------------------------------------
df_centered_normalized = read.csv("ALL_variance_partition_summary_centered_normalization.csv")
df_centered_normalized_1.0 = df_centered_normalized %>%
  filter(subset_value == 1.0) %>% 
  mutate(
    Prop_VA = Total_V_A / Total_V_G,
    Prop_VD = Total_V_D / Total_V_G,
    Prop_Cov = Total_cov / Total_V_G,
    
    Prop_ROH_VG = ROH_V_G / Total_V_G,
    Prop_NonROH_VG = Non_ROH_V_G / Total_V_G,
    Prop_ROH_Cov = 2*ROH_cov / Total_V_G
  ) %>%
  group_by(hXX, pop, tau, subset_value) %>%
  summarise(across(-simID, ~mean(.x, na.rm = TRUE)), .groups = 'drop') %>%
  select(-Total_check, -ROH_check) %>%
  mutate(
    hXX = case_when(
      hXX == "h00" ~ "Recessive",
      hXX == "h05" ~ "Additive",
      hXX == "h10" ~ "Dominant",
      hXX == 'h15' ~ "Mixed",
    ),
    hXX = factor(hXX, levels = c("Recessive", "Additive", "Dominant", "Mixed"))
  )%>%
  mutate(
    pop = case_when(
      pop == "p1" ~ "African",
      pop == "p2" ~ "European",
      pop == "p3" ~ "East Asian",
    ),
    pop = factor(pop, levels = c("African", "European", "East Asian"))
  )

df_centered_normalized_1.0_long = df_centered_normalized_1.0 %>%
  select(hXX, pop, tau, Prop_ROH_VG, Prop_NonROH_VG, Prop_ROH_Cov) %>%
  pivot_longer(
    cols = c(Prop_ROH_VG, Prop_NonROH_VG, Prop_ROH_Cov),
    names_to = "region",
    values_to = "proportion"
  )%>%
  mutate(region = factor(region, levels = c("Prop_ROH_VG", "Prop_NonROH_VG", "Prop_ROH_Cov")))
plot_centered_normalized_1.0 = ggplot(df_centered_normalized_1.0_long, aes(x = factor(tau), y = proportion, fill = region)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(hXX ~ pop) +
  scale_fill_manual(values = c(
    "Prop_ROH_VG" = "#A9DEF9",      
    "Prop_NonROH_VG" = "#90A955",   
    "Prop_ROH_Cov" = "#5E6472"     
  ),
  labels = c(
    "Prop_ROH_VG" = "ROH",
    "Prop_NonROH_VG" = "Non-ROH",
    "Prop_ROH_Cov" = "Covariance"
  )) +
  labs(
    x = "Tau",
    y = expression(bold("Proportion of V"[G]~"explained")),
    fill = "Region",
    title = "x = 1.0"
  ) +
  coord_cartesian(ylim = c(-0.05, 1.05)) +
  theme_bw() +
  theme(
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)) 
plot_centered_normalized_1.0
ggsave("Genetic_Variance_explained_by_ROH_centered_normalized_1.0.pdf", plot = plot_centered_normalized_1.0, width = 8, height = 8, device = cairo_pdf)

# # x=0.3(uncentered and normalized) ------------------------------------------
# df_uncentered_normalized = read.csv("ALL_variance_partition_summary_uncentered_normalization.csv")
# df_uncentered_normalized_0.3 = df_uncentered_normalized %>%
#   filter(subset_value == 0.3) %>% 
#   mutate(
#     Prop_VA = Total_V_A / Total_V_G,
#     Prop_VD = Total_V_D / Total_V_G,
#     Prop_Cov = Total_cov / Total_V_G,
#     
#     Prop_ROH_VG = ROH_V_G / Total_V_G,
#     Prop_NonROH_VG = Non_ROH_V_G / Total_V_G,
#     Prop_ROH_Cov = 2*ROH_cov / Total_V_G
#   ) %>%
#   group_by(hXX, pop, tau, subset_value) %>%
#   summarise(across(-simID, ~mean(.x, na.rm = TRUE)), .groups = 'drop') %>%
#   select(-Total_check, -ROH_check) %>%
#   mutate(
#     hXX = case_when(
#       hXX == "h00" ~ "Recessive",
#       hXX == "h05" ~ "Additive",
#       hXX == "h10" ~ "Dominant",
#       hXX == 'h15' ~ "Mixed",
#     ),
#     hXX = factor(hXX, levels = c("Recessive", "Additive", "Dominant", "Mixed"))
#   )%>%
#   mutate(
#     pop = case_when(
#       pop == "p1" ~ "African",
#       pop == "p2" ~ "European",
#       pop == "p3" ~ "East Asian",
#     ),
#     pop = factor(pop, levels = c("African", "European", "East Asian"))
#   )
# 
# df_uncentered_normalized_0.3_long = df_uncentered_normalized_0.3 %>%
#   select(hXX, pop, tau, Prop_ROH_VG, Prop_NonROH_VG, Prop_ROH_Cov) %>%
#   pivot_longer(
#     cols = c(Prop_ROH_VG, Prop_NonROH_VG, Prop_ROH_Cov),
#     names_to = "region",
#     values_to = "proportion"
#   )
# plot_uncentered_normalized_0.3 = ggplot(df_uncentered_normalized_0.3_long, aes(x = factor(tau), y = proportion, fill = region)) +
#   geom_bar(stat = "identity", position = "stack") +
#   facet_grid(hXX ~ pop) +
#   scale_fill_manual(values = c(
#     "Prop_ROH_VG" = "#66b3e6",      # 深蓝
#     "Prop_NonROH_VG" = "#ffb3b3",   # 橙色
#     "Prop_ROH_Cov" = "#c1471f"      # 绿色
#   ),
#   labels = c(
#     "Prop_ROH_VG" = "ROH",
#     "Prop_NonROH_VG" = "Non-ROH",
#     "Prop_ROH_Cov" = "Covariance"
#   )) +
#   labs(
#     x = "Tau",
#     y = expression(bold("Proportion of V"[G]~"explained")),
#     fill = "Region",
#     title = "x = 0.3 (Uncentered and Normalized)"
#   ) +
#   coord_cartesian(ylim = c(-0.05, 1.05)) +
#   theme_bw() +
#   theme(
#     axis.title.x = element_text(face = "bold"),
#     axis.title.y = element_text(face = "bold"),
#     legend.title = element_text(face = "bold"),
#     plot.title = element_text(hjust = 0.5, face = "bold"),
#     strip.text = element_text(face = "bold"),
#     axis.text.x = element_text(angle = 45, hjust = 1)) 
# ggsave("Genetic_Variance_explained_by_ROH_uncentered_normalized_0.3.pdf", plot = plot_uncentered_normalized_0.3, width = 8, height = 8, device = cairo_pdf)
# # x=1.0(uncentered and normalized) ------------------------------------------
# df_uncentered_normalized = read.csv("ALL_variance_partition_summary_uncentered_normalization.csv")
# df_uncentered_normalized_1.0 = df_uncentered_normalized %>%
#   filter(subset_value == 1.0) %>% 
#   mutate(
#     Prop_VA = Total_V_A / Total_V_G,
#     Prop_VD = Total_V_D / Total_V_G,
#     Prop_Cov = Total_cov / Total_V_G,
#     
#     Prop_ROH_VG = ROH_V_G / Total_V_G,
#     Prop_NonROH_VG = Non_ROH_V_G / Total_V_G,
#     Prop_ROH_Cov = 2*ROH_cov / Total_V_G
#   ) %>%
#   group_by(hXX, pop, tau, subset_value) %>%
#   summarise(across(-simID, ~mean(.x, na.rm = TRUE)), .groups = 'drop') %>%
#   select(-Total_check, -ROH_check) %>%
#   mutate(
#     hXX = case_when(
#       hXX == "h00" ~ "Recessive",
#       hXX == "h05" ~ "Additive",
#       hXX == "h10" ~ "Dominant",
#       hXX == 'h15' ~ "Mixed",
#     ),
#     hXX = factor(hXX, levels = c("Recessive", "Additive", "Dominant", "Mixed"))
#   )%>%
#   mutate(
#     pop = case_when(
#       pop == "p1" ~ "African",
#       pop == "p2" ~ "European",
#       pop == "p3" ~ "East Asian",
#     ),
#     pop = factor(pop, levels = c("African", "European", "East Asian"))
#   )
# 
# df_uncentered_normalized_1.0_long = df_uncentered_normalized_1.0 %>%
#   select(hXX, pop, tau, Prop_ROH_VG, Prop_NonROH_VG, Prop_ROH_Cov) %>%
#   pivot_longer(
#     cols = c(Prop_ROH_VG, Prop_NonROH_VG, Prop_ROH_Cov),
#     names_to = "region",
#     values_to = "proportion"
#   )
# plot_uncentered_normalized_1.0 = ggplot(df_uncentered_normalized_1.0_long, aes(x = factor(tau), y = proportion, fill = region)) +
#   geom_bar(stat = "identity", position = "stack") +
#   facet_grid(hXX ~ pop) +
#   scale_fill_manual(values = c(
#     "Prop_ROH_VG" = "#66b3e6",      # 深蓝
#     "Prop_NonROH_VG" = "#ffb3b3",   # 橙色
#     "Prop_ROH_Cov" = "#c1471f"      # 绿色
#   ),
#   labels = c(
#     "Prop_ROH_VG" = "ROH",
#     "Prop_NonROH_VG" = "Non-ROH",
#     "Prop_ROH_Cov" = "Covariance"
#   )) +
#   labs(
#     x = "Tau",
#     y = expression(bold("Proportion of V"[G]~"explained")),
#     fill = "Region",
#     title = "x = 1.0 (Uncentered and Normalized)"
#   ) +
#   coord_cartesian(ylim = c(-0.05, 1.05)) +
#   theme_bw() +
#   theme(
#     axis.title.x = element_text(face = "bold"),
#     axis.title.y = element_text(face = "bold"),
#     legend.title = element_text(face = "bold"),
#     plot.title = element_text(hjust = 0.5, face = "bold"),
#     strip.text = element_text(face = "bold"),
#     axis.text.x = element_text(angle = 45, hjust = 1)) 
# ggsave("Genetic_Variance_explained_by_ROH_uncentered_normalized_1.0.pdf", plot = plot_uncentered_normalized_1.0, width = 8, height = 8, device = cairo_pdf)
