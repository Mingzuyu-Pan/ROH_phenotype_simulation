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

hXX_values = c('h00', 'h05', 'h10', 'h15') 

cols_to_log = c("mean_A_frac", "mean_B_frac", "mean_C_frac", "mean_NONE_frac")

final_combined_summary = NULL 

for (hXX in hXX_values) {
  path = paste0(hXX, '_summary_frac_deleterious_normalized.xlsx')
  
  current_summary =read_excel(path)
  current_summary = current_summary[, -c(1, 5)]
  current_summary[cols_to_log] = log(current_summary[cols_to_log])
  current_summary$hXX = hXX
  final_combined_summary = rbind(final_combined_summary, current_summary)
}


final_combined_summary$A_change = final_combined_summary$mean_A_frac - final_combined_summary$mean_NONE_frac
final_combined_summary$B_change = final_combined_summary$mean_B_frac - final_combined_summary$mean_NONE_frac
final_combined_summary$C_change = final_combined_summary$mean_C_frac - final_combined_summary$mean_NONE_frac

final_combined_summary =  final_combined_summary %>% select(-mean_A_frac, -mean_B_frac, -mean_C_frac, -mean_NONE_frac)

stat_summary = final_combined_summary %>% group_by(hXX,pop,tau) %>%
  summarise(mean_change_A = mean(A_change),
            sd_change_A = sd(A_change),
            count_A = length(A_change),
            mean_change_B = mean(B_change),
            sd_change_B = sd(B_change),
            count_B = length(B_change),
            mean_change_C = mean(C_change),
            sd_change_C = sd(C_change),
            count_C = length(C_change),
            .groups = 'drop')
stat_summary$se_change_A = stat_summary$sd_change_A/sqrt(stat_summary$count_A)
stat_summary$se_change_B = stat_summary$sd_change_B/sqrt(stat_summary$count_B)
stat_summary$se_change_C = stat_summary$sd_change_C/sqrt(stat_summary$count_C)
stat_summary = stat_summary[,-c(5,6,8,9,11,12)]
stat_summary$tau =as.numeric(stat_summary$tau)/100

write.xlsx(
  x = stat_summary,
  file = "./Table3_change_in_per_unit_se_deleterious.xlsx"
)

stat_summary_long = stat_summary %>%
  pivot_longer(
    cols = contains("_change_"),
    names_to = "Metric_Category",
    values_to = "Value"
  ) %>%
  separate(
    col = Metric_Category,
    into = c("Metric", "ROH_type"),
    sep = "_change_",
    remove = TRUE,
    extra = "merge" 
  ) %>%
  mutate(
    Metric = gsub("_change", "", Metric) 
  )

id_cols = c("hXX", "pop", "tau", "ROH_type") 

stat_summary_plot = stat_summary_long %>%
  pivot_wider(
    id_cols = all_of(id_cols),
    names_from = Metric,
    values_from = Value
  )


pop_labels = c("p1" = "African","p2" = "European","p3" = "East Asian")
hXX_labels = c("h00" = "recessive", "h05" = "additive", "h10" = "dominant", "h15" = "mixed")
ROH_labels = c("A" = "Short ROH", "B" = "Medium ROH", "C" = "Long ROH")

p_test = ggplot(stat_summary_plot, aes(x = tau, y = mean, color = hXX, linetype = pop)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se, color = hXX), width = 0.04, alpha = 0.8,linetype = "solid") +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1) +
  facet_wrap(~ ROH_type,
             labeller = labeller(ROH_type = ROH_labels)
  ) +
  scale_color_manual(
    name = "Genetic Model",  
    values = c("h00" ='#B5134E', "h05" = '#1E88E5', "h10" = '#FFC107', "h15" = '#5AA99B'),
    labels = hXX_labels
  ) +
  scale_linetype_manual(
    name = "Population History",
    values = c("p1" = "solid", "p2" = "dashed", "p3" = "dotted"),
    labels = pop_labels) +
  theme_minimal()+
  theme(
    panel.grid.minor = element_blank(),
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
    breaks = seq(-1, 6, by = 1)
  )+
  labs(
    x = "Tau",  
    y = "Difference in per cM contribution between ROH and Non-ROH Regions",  
  )

print(p_test)
ggsave(
  filename = "change_in_fraction_deleterious_with_se.pdf",
  plot     = p_test,
  device   = cairo_pdf,
  width    = 6,
  height   = 8,
  units    = "in"
)

