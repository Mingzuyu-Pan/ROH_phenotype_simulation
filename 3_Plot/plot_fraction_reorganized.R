library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(showtext)
library(svglite)
library(ggpubr)
library(rstatix)
library(dplyr)
library(stringr)
library(openxlsx)

hXX_list = c('h00','h05','h10','h15')

all_means_list = list()
summary_log_list = list()
stat_list = list()
for (hXX in hXX_list) {
  if (hXX == 'h15'){
    name = 'mixed'
  } else if (hXX == 'h25'){
    name = 'mixed(ratio is based on chr1)'
  } else if (hXX == 'h00'){
    name = 'recessive'
  } else if (hXX == 'h05'){
    name = 'additive'
  } else if (hXX == 'h10'){
    name = 'dominant'
  }
path = paste0(hXX,'_summary_frac_deleterious_normalized.xlsx')
summary = read_excel(path)
summary = summary[,-1]
summary = summary[,-4]


cols_to_log = c("mean_A_frac", "mean_B_frac","mean_C_frac", "mean_NONE_frac")

# to get the log value
summary_new = summary
summary_new[cols_to_log] = log(summary[cols_to_log])
summary = summary_new

summary_new$hXX = hXX
summary_log_list[[hXX]] = summary_new

summary = summary %>%
  mutate(
    tau = as.numeric(as.character(tau)),  
  )
summary$tau = summary$tau / 100

tau_values = unique(summary$tau)
my_labels = c("p1" = "African","p2" = "European","p3" = "East Asian")

for (t in tau_values) {
  
summary_long = summary %>%pivot_longer(cols = starts_with("mean"), names_to = "metric", values_to = "value") 
summary_long = summary_long %>%filter(tau %in% t)

stat_test = summary_long %>%
  group_by(pop, tau) %>%
  #wilcox test and holm adjustment
  wilcox_test(value ~ metric, paired = TRUE, p.adjust.method = "holm") %>%
  add_significance(p.col = "p.adj", cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")) %>%
  add_xy_position(x = "pop", group = "metric", dodge = 0.7) %>%
  ungroup() %>%
  mutate(
    # to define the position of markers
    y.position = case_when(
      group2 == "mean_NONE_frac" & group1 == "mean_A_frac" ~ 1.45,
      group2 == "mean_NONE_frac" & group1 == "mean_B_frac" ~ 0.85,
      group2 == "mean_NONE_frac" & group1 == "mean_C_frac" ~ 0.25,
      
      group1 == "mean_A_frac"  & group2 == "mean_B_frac" ~ -10.50,
      group1 == "mean_A_frac"  & group2 == "mean_C_frac" ~ -11.75,
      group1 == "mean_B_frac"  & group2 == "mean_C_frac" ~ -11.15,
      
      TRUE ~ 3.0 
    )
  )


final_table <- stat_test %>%
  select(1, 2, 9, 10, 11, 13)

temp_name = paste0(hXX,'_',t)

stat_list[[temp_name]] <- final_table

plot_box = ggplot(summary_long, aes(x = pop, y = value, color = metric, fill = metric)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5, outlier.alpha = 0.5, position = position_dodge(width = 0.7), width = 0.5, size=0.1) +
  labs(
    title = name,
    x = "Tau",
    y = "Phenotype score per cM(log)",
    color = "Metric",
  ) +
  scale_color_manual(values = c("mean_A_frac" = "#fec590", "mean_B_frac" = "#fc8d59","mean_C_frac"="#d73027","mean_NONE_frac"="#4575b4"),labels = c("mean_A_frac" = "Short_frac", "mean_B_frac" = "Medium", "mean_C_frac" = "Long", "mean_NONE_frac" = "None")) +
 scale_fill_manual(name = "ROH Category", values = c("mean_A_frac" = "#fec590", "mean_B_frac" = "#fc8d59","mean_C_frac"="#d73027","mean_NONE_frac"="#4575b4"),labels = c("mean_A_frac" = "Short ROH", "mean_B_frac" = "Medium ROH", "mean_C_frac" = "Long ROH", "mean_NONE_frac" = "None ROH")) +
  

  theme(
    plot.title = element_text(family = "Arial", face = "bold", size = 16, hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(family = "Arial", face = "plain", size = 12),
    legend.position = "none"
  )+
  scale_x_discrete(labels = my_labels) +
  scale_y_continuous(
    limits = c(-12.5, 2),
    expand = expansion(add = c(0, 0))
  )
plot_box = plot_box +
  stat_pvalue_manual(
    stat_test %>% filter(y.position > 0), 
    label = "p.adj.signif",
    tip.length = 0.01, # define how marker looks
    hide.ns = FALSE,
    vjust = 0.65,
    family = "Arial"
  )

plot_box = plot_box +
  stat_pvalue_manual(
    stat_test %>% filter(y.position < 0), 
    label = "p.adj.signif",
    tip.length = -0.01, 
    vjust = 2.0,        
    hide.ns = FALSE,
    family = "Arial"
  )
print(plot_box)

ggsave(
  filename = paste0(t,"_frac_deleterious_",hXX,".pdf"),
  plot     = plot_box,
  device   = cairo_pdf,  
  width    = 3.5,
  height   = 6,
  units    = "in"
)




}

# combined_summary_new_log = do.call(rbind, summary_log_list)
# rownames(combined_summary_new_log) = NULL 
# write.csv(combined_summary_new_log,"fraction_error_bar_deleterious.csv", row.names = FALSE)



summary_long = summary %>%pivot_longer(cols = starts_with("mean"), names_to = "metric", values_to = "value") 
mean_table = summary_long %>%
  group_by(pop, tau, metric) %>%
  summarise(mean = mean(value, na.rm = TRUE),
            .groups = 'drop')
mean_table$hXX = hXX
all_means_list[[hXX]] = mean_table
}

final_summary_long_deleterious = dplyr::bind_rows(all_means_list)

# this calculation is decided for other following plotting and analysis
final_summary_wide_deleterious = final_summary_long_deleterious %>%
  pivot_wider(
    names_from = metric,     
    values_from = mean 
  )
write.csv(final_summary_wide_deleterious, "deleterious_mean_all_hXX.csv", row.names = FALSE)

stat_all <- bind_rows(stat_list, .id = "id_name")

stat_processed <- stat_all %>%
  mutate(
    sheet_group = str_extract(id_name, "^[^_]+") 
  )

stat_processed = stat_processed[,-1]

stat_final_list <- split(stat_processed, stat_processed$sheet_group)

write.xlsx(stat_final_list, file = "Table1_deleterious_statistical_results.xlsx")







hXX_list = c('h00','h05','h10','h15')

all_means_list = list()
summary_log_list = list()
stat_list = list()
for (hXX in hXX_list) {
  if (hXX == 'h15'){
    name = 'mixed'
  } else if (hXX == 'h25'){
    name = 'mixed(ratio is based on chr1)'
  } else if (hXX == 'h00'){
    name = 'recessive'
  } else if (hXX == 'h05'){
    name = 'additive'
  } else if (hXX == 'h10'){
    name = 'dominant'
  }
  path = paste0(hXX,'_summary_frac_neutral_normalized.xlsx')
  summary = read_excel(path)
  summary = summary[,-1]

  cols_to_log = c("mean_A_frac", "mean_B_frac","mean_C_frac", "mean_NONE_frac")

  # to get the log value
  summary_new = summary
  summary_new[cols_to_log] = log(summary[cols_to_log])
  summary = summary_new

  summary_new$hXX = hXX
  summary_log_list[[hXX]] = summary_new


  summary = summary %>%
    mutate(
      tau = as.numeric(as.character(tau)),  
      rho = as.numeric(as.character(rho))   
    )

  tau_values = unique(summary$tau)
  my_labels = c("p1" = "African","p2" = "European","p3" = "East Asian")

  for (t in tau_values) {

    summary_long = summary %>%pivot_longer(cols = starts_with("mean"), names_to = "metric", values_to = "value")
    summary_long = summary_long %>%filter(tau %in% t)

    stat_test = summary_long %>%
      group_by(pop, tau) %>%
      wilcox_test(value ~ metric, paired = TRUE, p.adjust.method = "holm") %>%
      add_significance(p.col = "p.adj", cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")) %>%
      add_xy_position(x = "pop", group = "metric", dodge = 0.7) %>%
      ungroup() %>%
      mutate(
        y.position = case_when(
          group2 == "mean_NONE_frac" & group1 == "mean_A_frac" ~ 1.05,
          group2 == "mean_NONE_frac" & group1 == "mean_B_frac" ~ 0.65,
          group2 == "mean_NONE_frac" & group1 == "mean_C_frac" ~ 0.25,
          
          group1 == "mean_A_frac"  & group2 == "mean_B_frac" ~ -5,
          group1 == "mean_A_frac"  & group2 == "mean_C_frac" ~ -5.8,
          group1 == "mean_B_frac"  & group2 == "mean_C_frac" ~ -5.4,
          TRUE ~ 3.0
        )
      )
    
    final_table <- stat_test %>%
      select(1, 2, 9, 10, 11, 13)
    
    temp_name = paste0(hXX,'_',t)
    
    stat_list[[temp_name]] <- final_table

    plot_box = ggplot(summary_long, aes(x = pop, y = value, color = metric, fill = metric)) +
      geom_boxplot(alpha = 0.7, outlier.size = 0.5, outlier.alpha = 0.5, position = position_dodge(width = 0.7), width = 0.5, size=0.1) +
      labs(
        title = name,
        x = "Tau",
        y = "Phenotype score per cM(log)",
        color = "Metric",
      ) +
      scale_color_manual(values = c("mean_A_frac" = "#fec590", "mean_B_frac" = "#fc8d59","mean_C_frac"="#d73027","mean_NONE_frac"="#4575b4"),labels = c("mean_A_frac" = "Short_frac", "mean_B_frac" = "Medium", "mean_C_frac" = "Long", "mean_NONE_frac" = "None")) +
      scale_fill_manual(name = "ROH Category", values = c("mean_A_frac" = "#fec590", "mean_B_frac" = "#fc8d59","mean_C_frac"="#d73027","mean_NONE_frac"="#4575b4"),labels = c("mean_A_frac" = "Short ROH", "mean_B_frac" = "Medium ROH", "mean_C_frac" = "Long ROH", "mean_NONE_frac" = "None ROH")) +


      theme(
        plot.title = element_text(family = "Arial", face = "bold", size = 16, hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(family = "Arial", face = "plain", size = 12),
        legend.position = "none"
      )+
      scale_x_discrete(labels = my_labels) +
      scale_y_continuous(
        limits = c(-6.5, 1.5),
        expand = expansion(add = c(0, 0))
      )
    plot_box = plot_box +
      stat_pvalue_manual(
        stat_test %>% filter(y.position > 0),
        label = "p.adj.signif",
        tip.length = 0.01,
        vjust = 0.7,
        hide.ns = FALSE,
        family = "Arial"
      )
    
    plot_box = plot_box +
      stat_pvalue_manual(
        stat_test %>% filter(y.position < 0),
        label = "p.adj.signif",
        tip.length = -0.01,
        vjust = 2.0,    
        hide.ns = FALSE,
        family = "Arial"
      )
    print(plot_box)

    ggsave(
      filename = paste0(t,"_frac_neutral_",hXX,".pdf"),
      plot     = plot_box,
      device   = cairo_pdf,
      width    = 3.5,
      height   = 6,
      units    = "in"
    )
  }

  # combined_summary_new_log = do.call(rbind, summary_log_list)
  # rownames(combined_summary_new_log) = NULL
  # write.csv(combined_summary_new_log,"fraction_error_bar_neutral.csv", row.names = FALSE)

  summary_long = summary %>%pivot_longer(cols = starts_with("mean"), names_to = "metric", values_to = "value")
  mean_table = summary_long %>%
    group_by(pop, tau, metric) %>%
    summarise(mean = mean(value, na.rm = TRUE),
    .groups = 'drop')
  mean_table$hXX = hXX
  all_means_list[[hXX]] = mean_table
}

final_summary_long_neutral = dplyr::bind_rows(all_means_list)

# this calculation is decided for other following plotting and analysis
final_summary_wide_neutral = final_summary_long_neutral %>%
  pivot_wider(
    names_from = metric,      
    values_from = mean  
  )

write.csv(final_summary_wide_neutral, "neutral_mean_all_hXX.csv", row.names = FALSE)
stat_all <- bind_rows(stat_list, .id = "id_name")

stat_processed <- stat_all %>%
  mutate(
    sheet_group = str_extract(id_name, "^[^_]+") 
  )

stat_processed = stat_processed[,-1]

stat_final_list <- split(stat_processed, stat_processed$sheet_group)

write.xlsx(stat_final_list, file = "Table2_neutral_statistical_results.xlsx")




