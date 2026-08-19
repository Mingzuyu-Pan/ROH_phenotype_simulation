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

data_100MB_chr1 = clean_empirical_data("All_Populations_chr1_100MB_ROH_Summary.csv", total_length = 127.8)

sim_h00 = clean_simulation_data("h00_simulated_ROH_Summary.csv", total_length = 127.85)
sim_h05 = clean_simulation_data("h05_simulated_ROH_Summary.csv", total_length = 127.85)
sim_h10 = clean_simulation_data("h10_simulated_ROH_Summary.csv", total_length = 127.85)


data_100MB_chr1_processed = data_100MB_chr1 %>%
    group_by(superpop, class) %>%
    summarise(
      mean_nROH = mean(nROH),
      se_nROH   = sd(nROH) / sqrt(n()),
      .groups = "drop"
    )%>%
  rename(pop = superpop)%>%
  mutate(name = "Empirical: Chr1 (100Mb)")

sim_h00_processed = sim_h00 %>%
  group_by(sim, pop, class) %>%
  summarise(
    mean_nROH_temp = mean(nROH),
    .groups = "drop"
  )%>%
  group_by(pop, class) %>%
  summarise(
    mean_nROH = mean(mean_nROH_temp),
    se_nROH   = sd(mean_nROH_temp) / sqrt(n()),
    .groups = "drop"
  )%>%
  mutate(name = "Simulated: Chr 1 (100 Mb, Recessive)")

sim_h05_processed = sim_h05 %>%
  group_by(sim, pop, class) %>%
  summarise(
    mean_nROH_temp = mean(nROH),
    .groups = "drop"
  )%>%
  group_by(pop, class) %>%
  summarise(
    mean_nROH = mean(mean_nROH_temp),
    se_nROH   = sd(mean_nROH_temp) / sqrt(n()),
    .groups = "drop"
  )%>%
  mutate(name = "Simulated: Chr 1 (100 Mb, Additive)")

sim_h10_processed = sim_h10 %>%
  group_by(sim, pop, class) %>%
  summarise(
    mean_nROH_temp = mean(nROH),
    .groups = "drop"
  )%>%
  group_by(pop, class) %>%
  summarise(
    mean_nROH = mean(mean_nROH_temp),
    se_nROH   = sd(mean_nROH_temp) / sqrt(n()),
    .groups = "drop"
  )%>%
  mutate(name = "Simulated: Chr 1 (100 Mb, Dominant)")

data_plot = bind_rows(data_100MB_chr1_processed, sim_h00_processed, sim_h05_processed, sim_h10_processed)%>%
  mutate(class = recode(class,
                        "A"     = "Short ROH",
                        "B"     = "Medium ROH",
                        "C"     = "Long ROH",
                        "total" = "Total ROH"
  ),
  class = factor(class, levels = c("Short ROH", "Medium ROH", "Long ROH", "Total ROH")),
  pop   = factor(pop, levels = c("AFR", "EUR", "EAS")))

plot = ggplot(data_plot, aes(x = pop, y = mean_nROH, fill = name)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(
    aes(ymin = mean_nROH - se_nROH, ymax = mean_nROH + se_nROH),
    position = position_dodge(width = 0.8),
    width = 0.2
  ) +
  facet_wrap(~ class, scales = "free_y") +
  labs(x = "Population", y = "Mean nROH", fill = NULL) +
  scale_fill_manual(values = c(
    "Empirical: Chr1 (100Mb)"                = "#2A9D8F",
    "Simulated: Chr 1 (100 Mb, Recessive)"   = "#E9C46A",
    "Simulated: Chr 1 (100 Mb, Additive)"    = "#F4A261",
    "Simulated: Chr 1 (100 Mb, Dominant)"     = "#E76F51"
  ))+
  theme_bw()
plot
ggsave(filename = "nROH_amplified.pdf", 
         plot = plot, width = 6.5, height = 8, bg = "white")
