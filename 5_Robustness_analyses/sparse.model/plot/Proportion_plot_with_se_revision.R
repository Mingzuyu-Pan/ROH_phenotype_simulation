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

df_long <- df_long %>%
  mutate(
    Region_type = factor(Region_type, levels = c("Short ROH", "Medium ROH", "Long ROH", "Non-ROH")),
    pop = factor(pop, levels = c("African", "European", "East Asian")),
    hXX = factor(hXX, levels = c("recessive", "additive", "dominant", "mixed")),)

df_ready = df_long %>% group_by(hXX,pop,tau,Region_type) %>%
  summarise(mean_proportion_plot = mean(proportion_mean),
            se_proportion_plot = sd(proportion_mean)/sqrt(500),
            .groups = 'drop')
color_palette = c("recessive" ='#B5134E', "additive" = '#1E88E5', "dominant" = '#FFC107', "mixed" = '#5AA99B')
line_styles = c("African" = "solid", "European" = "dashed", "East Asian" = "dotted")

tag_df <- data.frame(
  Region_type = factor(c("Short ROH", "Medium ROH", "Long ROH", "Non-ROH"),
                       levels = c("Short ROH", "Medium ROH", "Long ROH", "Non-ROH")),
  label = c("A", "B", "C", "D")
)

p_line <- ggplot(df_ready, aes(x = tau, y = mean_proportion_plot, 
                              color = hXX, linetype = pop, 
                              group = interaction(pop, hXX))) +
  
  geom_errorbar(aes(ymin = mean_proportion_plot - se_proportion_plot,
                    ymax = mean_proportion_plot + se_proportion_plot),
                width = 0.03, linewidth = 0.5) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1) +
  
  facet_wrap(~ Region_type, scales = "free_y", ncol = 2) +
  
  scale_color_manual(name = "Genetic model", values = color_palette) +
  scale_linetype_manual(name = "Population", values = line_styles) +
  
  scale_x_continuous(breaks = seq(0.7, 1.3, by = 0.1)) +
  
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  
  geom_text(data = tag_df, aes(x = -Inf, y = Inf, label = label),
            hjust = 0, vjust = -0.5, size = 5, fontface = "bold", 
            color = "black", inherit.aes = FALSE) +
  coord_cartesian(clip = "off") +
  
  labs(title = bquote(bold("x = " * .(x))),
       x = "Tau", 
       y = "Proportion of genotypic score explained") +
  
  theme_bw() +
  theme(
    text = element_text(family = "Arial"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12, margin = margin(b = 0)),
    strip.text = element_text(face = "bold", size = 11),
    strip.background = element_blank(), 
    
    axis.title.x = element_text(face = "bold", size = 11, margin = margin(t = 10,b = -10)),
    axis.title.y = element_text(face = "bold", size = 11, margin = margin(r = 10)),
    axis.text = element_text(color = "black", size = 10),
    legend.spacing.y = unit(-0.3, "cm"), 
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 10),
    legend.key.width = unit(1, "cm")
  ) +
  
  guides(
    linetype = guide_legend(order = 1, nrow = 1, override.aes = list(color = "black")),
    color = guide_legend(order = 2, nrow = 1)
  )

p_line
ggsave(
  filename = paste0('proportion_genotypic_score_explained_by_roh_with_se.x.',x,'.pdf'),
  plot     = p_line,
  device   = cairo_pdf,
  width    = 6.5,
  height   = 8,
  units    = "in"
)
