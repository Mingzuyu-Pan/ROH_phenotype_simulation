library(dplyr)
library(tidyr)
library(ggplot2)
df = read.csv("All_merged_results_number_of_causal.csv")
df_summary = df %>%
  group_by(hXX, pop, subset_value) %>%
  summarise(
    mean_ROH_single = mean(prop_ROH_single, na.rm = TRUE),
    mean_Non_ROH_single = mean(prpp_Non_ROH_single, na.rm = TRUE),
    mean_ROH_homo = mean(prop_ROH_homo, na.rm = TRUE),
    mean_Non_ROH_homo = mean(prpp_Non_ROH_homo, na.rm = TRUE),
    .groups = 'drop' 
  )%>%
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
df_single_long = df_summary %>%
  select(pop, hXX, subset_value, mean_ROH_single, mean_Non_ROH_single) %>%
  pivot_longer(
    cols = c(mean_ROH_single, mean_Non_ROH_single),
    names_to = "Region", 
    values_to = "Proportion"
  ) %>%
  mutate(Region = ifelse(Region == "mean_ROH_single", "ROH", "Non-ROH"))
# single causal allele (x = 0.3) --------------------------------------------
df_single_ready_0.3 =  df_single_long %>% filter(subset_value == 0.3) 
p_single_0.3 = ggplot(df_single_ready_0.3, aes(x = pop, y = Proportion, fill = Region)) +
  geom_col(width = 0.6) + 
  facet_wrap(~ hXX) +    
  scale_y_continuous() + 
  scale_fill_manual(values = c("ROH" = "#f2878a", "Non-ROH" = "#95b4cc")) + 
  theme_bw() +    
  labs(
    title = "x = 0.3",
    x = "Population",
    y = "Proportion of all causal variants",
    fill = "Genomic Region"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right",
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    legend.title = element_text(face = "bold") 
  )
p_single_0.3
ggsave("number_of_all_causal_0.3.pdf", plot = p_single_0.3, width = 6.5, height = 8.5, device = cairo_pdf)

# single causal allele (x = 1.0) --------------------------------------------
df_single_ready_1.0 =  df_single_long %>% filter(subset_value == 1.0) 
p_single_1.0 = ggplot(df_single_ready_1.0, aes(x = pop, y = Proportion, fill = Region)) +
  geom_col(width = 0.6) + 
  facet_wrap(~ hXX) +    
  scale_y_continuous() + 
  scale_fill_manual(values = c("ROH" = "#f2878a", "Non-ROH" = "#95b4cc")) + 
  theme_bw() +    
  labs(
    title = "x = 1.0",
    x = "Population",
    y = "Proportion of all causal variants",
    fill = "Genomic Region"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right",
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"), 
    legend.title = element_text(face = "bold") 
  )
p_single_1.0
ggsave("number_of_all_causal_1.0.pdf", plot = p_single_1.0, width = 6.5, height = 8.5, device = cairo_pdf)

# homo allele -----------------------------------------------------------
df_homo_long = df_summary %>%
  select(pop, hXX, subset_value, mean_ROH_homo, mean_Non_ROH_homo) %>%
  pivot_longer(
    cols = c(mean_ROH_homo, mean_Non_ROH_homo),
    names_to = "Region", 
    values_to = "Proportion"
  ) %>%
  mutate(Region = ifelse(Region == "mean_ROH_homo", "ROH", "Non-ROH"))

# homo causal allele (x = 0.3) --------------------------------------------
df_homo_ready_0.3 =  df_homo_long %>% filter(subset_value == 0.3) 
p_homo_0.3 = ggplot(df_homo_ready_0.3, aes(x = pop, y = Proportion, fill = Region)) +
  geom_col(width = 0.6) + 
  facet_wrap(~ hXX) +    
  scale_y_continuous() + 
  scale_fill_manual(values = c("ROH" = "#f2878a", "Non-ROH" = "#95b4cc")) + 
  theme_bw() +    
  labs(
    title = "x = 0.3",
    x = "Population",
    y = "Proportion of homozygous causal variants",
    fill = "Genomic Region"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right",
    axis.title.x = element_text(face = "bold"),  # X轴标题加粗
    axis.title.y = element_text(face = "bold"),  # Y轴标题加粗
    legend.title = element_text(face = "bold") 
  )
p_homo_0.3
ggsave("number_of_all_homo_0.3.pdf", plot = p_homo_0.3, width = 6.5, height = 8.5, device = cairo_pdf)

# single causal allele (x = 1.0) --------------------------------------------
df_homo_ready_1.0 =  df_homo_long %>% filter(subset_value == 1.0) 
p_homo_1.0 = ggplot(df_homo_ready_1.0, aes(x = pop, y = Proportion, fill = Region)) +
  geom_col(width = 0.6) + 
  facet_wrap(~ hXX) +    
  scale_y_continuous() + 
  scale_fill_manual(values = c("ROH" = "#f2878a", "Non-ROH" = "#95b4cc")) + 
  theme_bw(base_family = "Arial") +    
  labs(
    title = "x = 1.0",
    x = "Population",
    y = "Proportion of homozygous causal variants",
    fill = "Genomic Region"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right",
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    legend.title = element_text(face = "bold") 
  )
p_homo_1.0
ggsave("number_of_all_homo_1.0.pdf", plot = p_homo_1.0, width = 6.5, height = 8.5, device = cairo_pdf)
