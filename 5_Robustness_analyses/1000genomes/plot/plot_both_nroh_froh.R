library(ggplot2)
library(dplyr)
library(ggforce)
options(scipen = 999)


clean_empirical_data = function(file_path, total_length) {
  raw_data = read.csv(file_path)
  
  raw_data$superpop = factor(raw_data$superpop, levels = c("AFR", "EUR", "EAS"))
  raw_data_temp = raw_data %>%
    group_by(ind, pop, superpop) %>% 
    summarise(
      across(c(nROH, sROH_cM, sROH_bp), \(x) sum(x, na.rm = TRUE)), 
      .groups = "drop" 
    ) %>%
    mutate(class = "total")
  raw_data_final = bind_rows(raw_data, raw_data_temp) %>%
    mutate(class = factor(class, levels = c("A", "B", "C", "total"))) %>%
    arrange(ind, class)
  raw_data_final$FROH_cM = raw_data_final$sROH_cM/total_length
  return(raw_data_final)
}

clean_simulation_data = function(file_path, total_length) {
  raw_data = read.csv(file_path)
  raw_data_clean = raw_data %>%
    mutate(pop = case_when(
      pop == "p1" ~ "AFR",
      pop == "p2" ~ "EUR",
      pop == "p3" ~ "EAS",
      TRUE ~ pop
    ))
  raw_data_clean$pop = factor(raw_data_clean$pop, levels = c("AFR", "EUR", "EAS"))
  raw_data_totals = raw_data_clean %>%
    group_by(sim, ind, pop) %>% 
    summarise(
      across(c(nROH, sROH_cM, sROH_bp), \(x) sum(x, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    mutate(class = "total")
  raw_data_final = bind_rows(raw_data_clean, raw_data_totals) %>%
    mutate(class = factor(class, levels = c("A", "B", "C", "total"))) %>%
    arrange(sim, ind, pop, class)
  raw_data_final$FROH_cM = raw_data_final$sROH_cM/total_length

  return(raw_data_final)
}
generate_summary_plots = function(data, pop_var, label) {

  save_plot = function(plot_obj, title) {
    clean_name <- gsub("^_|_$", "", gsub("[^[:alnum:]]+", "_", title))
    ggsave(filename = paste0(clean_name, ".pdf"), plot = plot_obj, width = 6, height = 4, bg = "white")
  }
  
  # Total sum of nROH
  title1 = paste0("Total sum of nROH - ", label)
  p1 = ggplot(data, aes(x = .data[[pop_var]], y = nROH, fill = .data[[pop_var]])) +
    stat_summary(fun = sum, geom = "bar", width = 0.3, color = "black", alpha = 0.8) + 
    facet_wrap(~ class, scales = "free_y") +         
    theme_bw() +                                   
    labs(title = title1, x = "Population", y = "Total Number of ROH (Sum of nROH)") +
    theme(legend.position = "none", strip.background = element_rect(fill = "grey90"), plot.title = element_text(hjust = 0.5))
  save_plot(p1, title1)
  
  # Average nROH per individual
  title2 = paste0("Average nROH per individual - ", label)
  p2 = ggplot(data, aes(x = .data[[pop_var]], y = nROH, fill = .data[[pop_var]])) +
    stat_summary(fun = mean, geom = "bar", width = 0.3, color = "black", alpha = 0.8) + 
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) + 
    facet_wrap(~ class, scales = "free_y") +         
    theme_bw() +                                   
    labs(title = title2, x = "Population", y = "Average Number of ROH per Individual (Mean nROH)") + 
    theme(legend.position = "none", strip.background = element_rect(fill = "grey90"), plot.title = element_text(hjust = 0.5))
  save_plot(p2, title2)
  
  # FROH (Boxplot)
  title3 = paste0("FROH - ", label)
  data_total = data %>% filter(class == "total")
  p3 = ggplot(data_total, aes(x = .data[[pop_var]], y = FROH_cM, fill = .data[[pop_var]])) +
    geom_boxplot(width = 0.5, outlier.size = 1, outlier.alpha = 0.3, outlier.stroke = 0, alpha = 0.8) + 
    scale_y_sqrt(limits = c(0, 1.0)) + 
    stat_summary(fun = median, geom = "text", aes(label = sprintf("%.3f", after_stat(y))), vjust = -1.2, size = 3.5, color = "red") + 
    theme_bw() +                                   
    labs(title = title3, x = "Population", y = "Total FROH") +
    theme(legend.position = "none", strip.background = element_rect(fill = "grey90"), plot.title = element_text(hjust = 0.5))
  save_plot(p3, title3)
}

data_all_chr = clean_empirical_data("All_Populations_ROH_Summary.csv", total_length = 3545.8255)
data_whole_chr1 = clean_empirical_data("All_Populations_chr1_ROH_Summary.csv", total_length = 286.2792)
data_100MB_chr1 = clean_empirical_data("All_Populations_chr1_100MB_ROH_Summary.csv", total_length = 127.8)

generate_summary_plots(data_all_chr,    pop_var = "superpop", label = "all chr + empirical")
generate_summary_plots(data_whole_chr1, pop_var = "superpop", label = "whole chr1 + empirical")
generate_summary_plots(data_100MB_chr1, pop_var = "superpop", label = "100 Mbps chr1 + empirical")


sim_h00 = clean_simulation_data("h00_simulated_ROH_Summary.csv", total_length = 127.85)
sim_h05 = clean_simulation_data("h05_simulated_ROH_Summary.csv", total_length = 127.85)
sim_h10 = clean_simulation_data("h10_simulated_ROH_Summary.csv", total_length = 127.85)

generate_summary_plots(sim_h00, pop_var = "pop", label = "100 Mbps chr1 + simulation + recessive")
generate_summary_plots(sim_h05, pop_var = "pop", label = "100 Mbps chr1 + simulation + additive")
generate_summary_plots(sim_h10, pop_var = "pop", label = "100 Mbps chr1 + simulation + dominant")
