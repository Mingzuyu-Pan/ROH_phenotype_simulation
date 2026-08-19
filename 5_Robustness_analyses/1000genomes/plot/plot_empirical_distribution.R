library(ggplot2)
library(dplyr)
library(ggforce)
options(scipen = 999)


generate_empirical_distribution = function(file_path, label) {
  raw_data = read.csv(file_path)
  
  raw_data$superpop = factor(raw_data$superpop, levels = c("AFR", "EUR", "EAS"))
  raw_data = raw_data[raw_data$cM_length > 0, ]
  raw_data = raw_data[raw_data$bp_length > 0, ]

  
  build_plot <- function(data, x_var, x_label, title, x_limits) {
    p <- ggplot(data, aes(x = {{ x_var }}, fill = superpop)) +
      geom_histogram(aes(y = after_stat(count / ave(count, PANEL, FUN = sum))), 
                     binwidth = 0.1,
                     color = "white", alpha = 0.8) +
      scale_y_sqrt(limits = c(0, 0.15)) +
      scale_x_log10() + 
      coord_cartesian(xlim = x_limits) +
      facet_wrap(~ superpop, ncol = 1, scales = "free_y") + 
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
  
  build_plot(data = raw_data, 
             x_var = cM_length, 
             x_label = "ROH Length (cM) [Log10 Scale]", 
             title = paste0("ROH Length Spectrum by Population (cM) + ", label, " + empirical"),
             x_limits = c(0.00001, 12))
  
  build_plot(data = raw_data, 
             x_var = bp_length/1e6, 
             x_label = "ROH Length (Mb) [Log10 Scale]", 
             title = paste0("ROH Length Spectrum by Population (Mb) + ", label, " + empirical"),
             x_limits = c(0.00001, 12)) 
}

generate_empirical_distribution("All_Populations_ROH_Summary_distr.csv.gz", "all chr")
generate_empirical_distribution("All_Populations_chr1_ROH_Summary_distr.csv.gz", "whole chr1")
generate_empirical_distribution("All_Populations_chr1_100MB_ROH_Summary_distr.csv.gz","100 Mb chr1")
