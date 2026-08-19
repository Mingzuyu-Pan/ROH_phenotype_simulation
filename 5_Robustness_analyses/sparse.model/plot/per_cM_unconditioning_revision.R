library(tidyverse)
library(ggplot2) 
library(dplyr)
library(patchwork)
options(scipen = 999)
env_list = c('1.0')
x_list = c('1.0', '0.3')
env = '1.0'
x = '0.3'

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

my_colors <- c("Short ROH" = "#fec590",   
               "Medium ROH" = "#fc8d59",  
               "Long ROH" = "#d73027",    
               "Non-ROH" = "#4575b4")   

plot_panel_uncondition <- function(data, tau_val, show_x_title = FALSE, show_y_title = FALSE, 
                                 show_x_text = FALSE, tag_text = "", custom_margin = margin(5, 5, 5, 5)) {
  
  df_sub <- data %>% filter(tau == tau_val)
  
  stat_test <- df_sub %>%
    group_by(pop, hXX) %>%
    wilcox_test(per_cM_all_log ~ Region_type, paired = TRUE, p.adjust.method = "holm") %>%
    add_significance(p.col = "p.adj", cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                     symbols = c("***", "**", "*", "ns")) %>%
    add_xy_position(x = "pop", group = "Region_type", dodge = 0.8) %>%
    ungroup() %>%
    mutate(
      y.position = case_when(
        group2 == "Non-ROH" & group1 == "Short ROH"  ~ 0.7,
        group2 == "Non-ROH" & group1 == "Medium ROH" ~ 0,
        group2 == "Non-ROH" & group1 == "Long ROH"   ~ -0.7,
        group1 == "Short ROH"  & group2 == "Medium ROH" ~ -13.8,
        group1 == "Short ROH"  & group2 == "Long ROH"   ~ -14.7,
        group1 == "Medium ROH" & group2 == "Long ROH"   ~ -12.9,
        TRUE ~ 0.5
      )
    )
  
  x_title <- if(show_x_title) "Population history" else NULL
  y_title <- if(show_y_title) "Unconditioning phenotype score per cM (ln)" else NULL
  
  p <- ggplot(df_sub, aes(x = pop, y = per_cM_all_log, fill = Region_type)) +
    geom_boxplot(
      color = "grey7",    
      alpha = 0.7, outlier.size = 0.3, outlier.alpha = 0.3, position = position_dodge(width = 0.8), width = 0.7, size=0.1) +
    scale_fill_manual(values = my_colors) +
    facet_wrap(~ hXX, nrow = 1) +
    labs(
      tag = tag_text, 
      title = paste0("Tau = ", tau_val),
      x = x_title, 
      y = y_title,
      fill = "Region type"
    ) +
    scale_y_continuous(
      limits = c(-16, 2),
      expand = expansion(add = c(0, 0))
    )+
    theme_minimal() +
    theme(
      plot.tag = element_text(size = 18, face = "bold"),
      plot.tag.position = c(0.075, 0.93),
      
      plot.title = element_text(hjust = 0.5, face = "bold", size = 10, margin = margin(b = 0)),
      strip.text = element_text(face = "bold", size = 8),
      strip.background = element_blank(),
      axis.text.y = element_text(color = "black"),
      axis.title.x = element_text(face = "bold", size = 10, margin = margin(t = 10)),
      axis.title.y = element_text(face = "bold", size = 10),
      panel.background = element_rect(fill = "#F2F2F2", color = NA),
      panel.grid.major = element_line(color = "white", linewidth = 0.6),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "grey70", fill = NA, linewidth = 0.5),
      panel.spacing = unit(0, "cm"),
      text = element_text(family = "Arial"),
      plot.margin = custom_margin
    )
  
  p <- p +
    stat_pvalue_manual(
      stat_test %>% filter(y.position > -1),
      label = "p.adj.signif",
      tip.length = 0.01,
      hide.ns = FALSE,
      vjust = 0.65,
      family = "Arial",
      size = 2.5
    )
  
  p <- p +
    stat_pvalue_manual(
      stat_test %>% filter(y.position < -1),
      label = "p.adj.signif",
      tip.length = -0.01,
      vjust = 1.65,
      hide.ns = FALSE,
      family = "Arial",
      size = 2.5
    )
  
  if (show_x_text) {
    p <- p + theme(
      axis.ticks.x = element_blank(),
      axis.ticks.length.x = unit(0, "pt"),
      
      axis.text.x = element_text(
        angle = 45, 
        hjust = 1, 
        vjust = 1,              
        color = "black",
        margin = margin(t = -25)
      )
    )
  } else {
    p <- p + theme(
      axis.text.x = element_blank(), 
      axis.ticks.x = element_blank(),
      axis.ticks.length.x = unit(0, "pt")
    )
  }
  
  return(p)
}

pD <- plot_panel_uncondition(df_long, tau_val = 0.7, 
                           show_x_title = FALSE, show_y_title = FALSE, show_x_text = FALSE, tag_text = "A",
                           custom_margin = margin(t = 5, r = 10, b = 3, l = 15))

pE <- plot_panel_uncondition(df_long, tau_val = 1.0, 
                           show_x_title = FALSE, show_y_title = TRUE,  show_x_text = FALSE, tag_text = "B",
                           custom_margin = margin(t = 5, r = 10, b = 3, l = 15))

pF <- plot_panel_uncondition(df_long, tau_val = 1.3, 
                           show_x_title = TRUE,  show_y_title = FALSE, show_x_text = TRUE, tag_text = "C",
                           custom_margin = margin(t = 5, r = 10, b = 3, l = 15))

final_plot_uncondition <- pD / pE / pF + 
  plot_layout(guides = "collect")+plot_annotation(
    title = bquote(bold("x = " * .(x))),
    theme = theme(
      plot.title = element_text(hjust = 0.55, size = 12, margin = margin(b = 0)),
      plot.margin = margin(t = 0.75, r =0.25, b = 0.25, l = 0.25)
    )
  )  & 
  theme(
    legend.position = "bottom", 
    legend.title = element_text(face = "bold", size = 10),
    legend.box.margin = margin(t = -10, r = 0, b = 0, l = 0),

  )

final_plot_uncondition
ggsave(
  filename = paste0('per_cM_unconditioning.x.',x,'.pdf'),
  plot     = final_plot_uncondition,
  device   = cairo_pdf,
  width    = 6.5,
  height   = 8,
  units    = "in"
)

all_stat <- df_long %>%
  group_by(pop, tau, hXX) %>%
  wilcox_test(per_cM_all_log ~ Region_type, paired = TRUE, p.adjust.method = "holm") %>%
  add_significance(p.col = "p.adj", cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                   symbols = c("***", "**", "*", "ns"))

write.csv(all_stat, paste0('statistical_result_per_cM_unconditioning.x.',x,'.csv'), row.names = FALSE)
