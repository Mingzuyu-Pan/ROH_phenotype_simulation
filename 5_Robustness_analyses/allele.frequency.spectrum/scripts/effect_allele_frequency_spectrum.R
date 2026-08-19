library(tidyverse)
library(scales) 
hXX_list = c('h00','h05','h10','h15')

for (hXX in hXX_list) {
  if (hXX == 'h00'){
    name = 'Recessive'
  } else if (hXX == 'h05'){
    name = 'Additive'
  } else if (hXX == 'h10'){
    name = 'Dominant'
  }else if (hXX == 'h15'){
    name = 'Mixed'
  }
df = read.csv(paste0('effect_allele_frequency_spectrum_casual_effect_',hXX,'.subset.1.0.csv'))
df_long = df %>%
  pivot_longer(
    cols = -freq_bin, 
    names_to = c("effect", "pop", ".value"), 
    names_pattern = "effect_(.*)_(p[123])_(.*)"
  ) %>%
  filter(!is.na(mean)) 

df_long = df_long %>%
  mutate(
    left_edge_str = str_extract(freq_bin, "[0-9e.-]+"),
    left_edge_num = as.numeric(left_edge_str),
    log10_freq = log10(left_edge_num) 
  )%>%
  mutate(
    pop = case_when(
      pop == "p1" ~ "African",
      pop == "p2" ~ "European",
      pop == "p3" ~ "East Asian"
    ),
    pop = factor(pop, levels = c("African", "European", "East Asian"))
  )%>%
  mutate(
    Tau_label = factor(paste0("Tau = ", effect), 
                       levels = c("Tau = 0.7", "Tau = 1.0", "Tau = 1.3"))
  )


plot = ggplot(df_long, aes(x = log10_freq, y = mean, fill = pop)) +
  geom_col(width = 0.118, alpha = 0.8) +
  geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem), 
                width = 0.03, alpha = 0.7, color = "black") +
  facet_grid(Tau_label ~ pop) +
  scale_fill_manual(
    values = c(
      "African" = "#fba0a8", 
      "European"   = "#f8ce4f",  
      "East Asian" = "#90a9f6"   
    ),
  ) +
  ylim(0, 0.065) +
  
  scale_x_continuous(
    breaks = -6:0, 
    labels = label_math(10^.x) 
  ) +
  
  theme_bw(base_size = 12.5, base_family = "Arial") +
  labs(
    title = paste0(name, ' (x = 1.0)'),
    x = "Allele Frequency (log10 scale)",
    y = "Mean Effect",
    fill = "Population"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    strip.text = element_text(size = 12, face = "bold"),
    legend.title = element_text(face = "bold", size = 12),
    axis.text.x = element_text(size = 11), 
    axis.text.y = element_text(size = 10),
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )
print(plot)
ggsave(paste0('effect_allele_frequency_spectrum_',hXX,'.subset.1.0.pdf'), plot = plot, width = 9, height = 8, device = cairo_pdf)
}

for (hXX in hXX_list) {
  if (hXX == 'h00'){
    name = 'Recessive'
  } else if (hXX == 'h05'){
    name = 'Additive'
  } else if (hXX == 'h10'){
    name = 'Dominant'
  }else if (hXX == 'h15'){
    name = 'Mixed'
  }
  df = read.csv(paste0('effect_allele_frequency_spectrum_casual_effect_',hXX,'.subset.0.3.csv'))
  df_long = df %>%
    pivot_longer(
      cols = -freq_bin, 
      names_to = c("effect", "pop", ".value"), 
      names_pattern = "effect_(.*)_(p[123])_(.*)"
    ) %>%
    filter(!is.na(mean)) 
  
  df_long = df_long %>%
    mutate(
      left_edge_str = str_extract(freq_bin, "[0-9e.-]+"),
      left_edge_num = as.numeric(left_edge_str),
      log10_freq = log10(left_edge_num) 
    )%>%
    mutate(
      pop = case_when(
        pop == "p1" ~ "African",
        pop == "p2" ~ "European",
        pop == "p3" ~ "East Asian"
      ),
      pop = factor(pop, levels = c("African", "European", "East Asian"))
    )%>%
    mutate(
      Tau_label = factor(paste0("Tau = ", effect), 
                         levels = c("Tau = 0.7", "Tau = 1.0", "Tau = 1.3"))
    )
  
  
  plot = ggplot(df_long, aes(x = log10_freq, y = mean, fill = pop)) +
    geom_col(width = 0.118, alpha = 0.8) +
    geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem), 
                  width = 0.03, alpha = 0.7, color = "black") +
    facet_grid(Tau_label ~ pop) +
    scale_fill_manual(
      values = c(
        "African" = "#fba0a8", 
        "European"   = "#f8ce4f",  
        "East Asian" = "#90a9f6"   
      ),
    ) +
    ylim(0, 0.065) +
    
    scale_x_continuous(
      breaks = -6:0, 
      labels = label_math(10^.x) 
    ) +
    
    theme_bw(base_size = 12.5, base_family = "Arial") +
    labs(
      title = paste0(name, ' (x = 0.3)'),
      x = "Allele Frequency (log10 scale)",
      y = "Mean Effect",
      fill = "Population"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      strip.text = element_text(size = 12, face = "bold"),
      legend.title = element_text(face = "bold", size = 12),
      axis.text.x = element_text(size = 11), 
      axis.text.y = element_text(size = 10),
      axis.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
  print(plot)
  ggsave(paste0('effect_allele_frequency_spectrum_',hXX,'.subset.0.3.pdf'), plot = plot, width = 9, height = 8, device = cairo_pdf)
}

