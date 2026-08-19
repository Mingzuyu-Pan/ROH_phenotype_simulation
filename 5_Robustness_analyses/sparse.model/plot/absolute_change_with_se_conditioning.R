library(tidyverse)
library(ggplot2) 
library(dplyr)
library(patchwork)
options(scipen = 999)
env_list = c('1.0')
x_list = c('1.0', '0.3')
env = '1.0'
x = '1.0'

df_raw = read_csv(paste0('all_simulations_roh_pheno_vcf_revision.subset.',x,'.csv'))
df_long <- df_raw %>%
  pivot_longer(
    cols = -c(pop, tau, hXX, simID),
    names_to = c("Region_type", ".value"),
    names_pattern = "^(A|B|C|Non-ROH)_(.*)$"
  ) %>%
  mutate(
    Region_type = case_when(
      Region_type == "A" ~ "Short ROH",
      Region_type == "B" ~ "Medium ROH",
      Region_type == "C" ~ "Long ROH",
      Region_type == "Non-ROH" ~ "Non-ROH",
      TRUE ~ Region_type
    ),
    per_cM_all_log = log(per_cM_mean_all),
    per_cM_nonzero_log = log(per_cM_mean_nonzero),
    hXX = case_when(
      hXX == "h00" ~ "recessive",
      hXX == "h05" ~ "additive",
      hXX == "h10" ~ "dominant",
      hXX == "h15" ~ "mixed",
      TRUE ~ hXX 
    ),
    pop = case_when(
      pop == "p1" ~ "African",
      pop == "p2" ~ "European",
      pop == "p3" ~ "East Asian",
      TRUE ~ pop
    ),
    Region_type = factor(Region_type, levels = c("Short ROH", "Medium ROH", "Long ROH", "Non-ROH")),
    pop = factor(pop, levels = c("African", "European", "East Asian")),
    hXX = factor(hXX, levels = c("recessive", "additive", "dominant", "mixed")),)

df_diff <- df_long %>%
  group_by(pop, tau, hXX, simID) %>%
  mutate(
    log_per_cM_all = log(per_cM_mean_all),
    log_per_cM_nonzero = log(per_cM_mean_nonzero),
    diff_log_all = log_per_cM_all - log_per_cM_all[Region_type == "Non-ROH"],
    diff_log_nonzero = log_per_cM_nonzero - log_per_cM_nonzero[Region_type == "Non-ROH"]
  ) %>%
  ungroup() 

df_diff_filtered = df_diff[df_diff$Region_type != "Non-ROH", ]

df_ready = df_diff_filtered %>% group_by(hXX,pop,tau,Region_type) %>%
  summarise(mean_diff_log_nonzero = mean(diff_log_nonzero),
            se_diff_log_nonzero = sd(diff_log_nonzero)/sqrt(500),
            .groups = 'drop')

p_test = ggplot(df_ready, aes(x = tau, y = mean_diff_log_nonzero, color = hXX, linetype = pop)) +
  geom_errorbar(aes(ymin = mean_diff_log_nonzero - se_diff_log_nonzero, ymax =  mean_diff_log_nonzero + se_diff_log_nonzero, color = hXX), width = 0.04, alpha = 0.8,linetype = "solid") +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1) +
  facet_wrap(~ Region_type) +
  scale_color_manual(
    name = "Genetic Model",  
    values = c("recessive" ='#B5134E', "additive" = '#1E88E5', "dominant" = '#FFC107', "mixed" = '#5AA99B')
  ) +
  scale_linetype_manual(
    name = "Population History",
    values = c("African" = "solid", "European" = "dashed", "East Asian" = "dotted")) +
  theme_minimal()+
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12, margin = margin(b = 0)),
    text = element_text(family = "Arial"),
    axis.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    legend.title = element_text(
      face = "bold",
      size = 10,
      color = "black" 
    )
  )+
  scale_x_continuous(
    breaks = seq(0.7, 1.3, by = 0.1)
  )+
  scale_y_continuous(
    limits = c(-3.5,6.5),
    breaks = seq(-3.5, 6.5, by = 1)
  )+
  labs(
    title = bquote(bold("x = " * .(x))),
    x = "Tau",  
    y = "Difference in per cM contribution between ROH and Non-ROH Regions (Conditional)",  
  )

print(p_test)
ggsave(
  filename = paste0('change_per_cM_conditional.x.',x,'.pdf'),
  plot     = p_test,
  device   = cairo_pdf,
  width    = 6.5,
  height   = 8,
  units    = "in"
)
