library(arrow)
library(dplyr)
library(ggplot2)
options(scipen = 999)

generate_roh_plot <- function(file_path, sim_model) {
  dataset = open_dataset(file_path)
  plot_data = dataset %>%
    select(pop, cM_length, bp_length) %>% 
    collect() %>%
    mutate(pop = case_when(
      pop == "p1" ~ "AFR",
      pop == "p2" ~ "EUR",
      pop == "p3" ~ "EAS",
      TRUE ~ pop
    )) %>%
    filter(cM_length > 0, bp_length > 0)
  
  plot_data$pop = factor(plot_data$pop, levels = c("AFR", "EUR", "EAS"))
  
  build_plot <- function(data, x_var, x_label, title, x_limits) {
    p <- ggplot(data, aes(x = {{ x_var }}, fill = pop)) +
      geom_histogram(aes(y = after_stat(count / ave(count, PANEL, FUN = sum))), 
                     binwidth = 0.1,
                     color = "white", alpha = 0.8) +
      #scale_y_sqrt(limits = c(0, 0.15)) +
      scale_y_log10() +
      scale_x_log10() + 
      #coord_cartesian(xlim = x_limits) +
      facet_wrap(~ pop, ncol = 1, scales = "free_y") + 
      theme_bw() +
      labs(title = title,
           x = x_label,
           y = "Frequency (Proportion)") +
      theme(legend.position = "none",
            strip.background = element_rect(fill = "grey90"), 
            plot.title = element_text(hjust = 0.5))
    
    raw_name <- gsub("[^[:alnum:]]+", "_", title)
    clean_name <- gsub("^_|_$", "", raw_name) 
    
    ggsave(filename = paste0(clean_name, ".pdf"), 
           plot = p, width = 6, height = 4, bg = "white")
  }
  
  build_plot(data = plot_data, 
             x_var = cM_length, 
             x_label = "ROH Length (cM) [Log10 Scale]", 
             title = paste0("ROH Length Spectrum by Population (cM) - 100MB chr1 + simulation + ", sim_model),
             x_limits = c(0.00001, 12))
  
  build_plot(data = plot_data, 
             x_var = bp_length/1e6, 
             x_label = "ROH Length (Mb) [Log10 Scale]", 
             title = paste0("ROH Length Spectrum by Population (Mb) - 100MB chr1 + simulation + ", sim_model),
             x_limits = c(0.00001, 12)) 
}

generate_roh_plot("./h00_simulated_ROH_Summary.parquet_dataset", "recessive")
generate_roh_plot("./h05_simulated_ROH_Summary.parquet_dataset", "additive")
generate_roh_plot("./h10_simulated_ROH_Summary.parquet_dataset", "dominant")
