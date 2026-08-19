library(arrow)
library(dplyr)
library(ggplot2)
library(scales) 
options(scipen = 999)

load_simulation_data <- function(file_path, sim_model) {
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
  plot_data$condition = sim_model
  plot_data$pop_plotting = plot_data$pop
  plot_data = plot_data %>% 
    select(cM_length, bp_length, condition, pop_plotting)
  return(plot_data)
}

load_empirical_data = function(file_path, label) {
  raw_data = read.csv(file_path)
  
  raw_data$superpop = factor(raw_data$superpop, levels = c("AFR", "EUR", "EAS"))
  raw_data = raw_data[raw_data$cM_length > 0, ]
  raw_data = raw_data[raw_data$bp_length > 0, ]
  raw_data$condition = label
  raw_data$pop_plotting = raw_data$superpop
  plot_data = raw_data %>% 
    select(cM_length, bp_length, condition, pop_plotting)
  return(plot_data)
}

extract_official_ccdf_cM = function(data) {
  
  data_AFR = data %>% filter(pop_plotting == "AFR")
  data_EUR = data %>% filter(pop_plotting == "EUR")
  data_EAS = data %>% filter(pop_plotting == "EAS")
  
  ccdf_cal = function(data){
    current_ecdf = ecdf(data$cM_length)
    x_values = knots(current_ecdf)
    cdf_values = current_ecdf(x_values)
    y_values = 1 - cdf_values
    result = data.frame(
      pop_plotting = data$pop_plotting[1],
      cM_length = x_values,
      ccdf_y = y_values,
      condition = data$condition[1])
    return(result)
  }
  data_AFR_processed = ccdf_cal(data_AFR)
  data_EUR_processed = ccdf_cal(data_EUR)
  data_EAS_processed = ccdf_cal(data_EAS)
  final = bind_rows(data_AFR_processed, data_EUR_processed, data_EAS_processed)
  # 把三个人群的坐标拼起来
  return(final)
}

extract_official_ccdf_Mb = function(data) {
  
  data_AFR = data %>% filter(pop_plotting == "AFR")
  data_EUR = data %>% filter(pop_plotting == "EUR")
  data_EAS = data %>% filter(pop_plotting == "EAS")
  
  ccdf_cal = function(data){
    current_ecdf = ecdf(data$bp_length/1e6)
    x_values = knots(current_ecdf)
    cdf_values = current_ecdf(x_values)
    y_values = 1 - cdf_values
    result = data.frame(
      pop_plotting = data$pop_plotting[1],
      Mb_length = x_values,
      ccdf_y = y_values,
      condition = data$condition[1])
    return(result)
  }
  data_AFR_processed = ccdf_cal(data_AFR)
  data_EUR_processed = ccdf_cal(data_EUR)
  data_EAS_processed = ccdf_cal(data_EAS)
  final = bind_rows(data_AFR_processed, data_EUR_processed, data_EAS_processed)
  # 把三个人群的坐标拼起来
  return(final)
}

empirical_all_autosomes = load_empirical_data("All_Populations_ROH_Summary_distr.csv.gz", "Empirical: All Autosomes")
empirical_whole_chr1 = load_empirical_data("All_Populations_chr1_ROH_Summary_distr.csv.gz", "Empirical: Whole Chr 1")
empirical_100Mb_chr1 = load_empirical_data("All_Populations_chr1_100MB_ROH_Summary_distr.csv.gz","Empirical: Chr1 (100 Mb)")
simulation_h00 = load_simulation_data("./h00_simulated_ROH_Summary.parquet_dataset", "Simulated: Chr 1 (100 Mb, Recessive)")
simulation_h05 = load_simulation_data("./h05_simulated_ROH_Summary.parquet_dataset", "Simulated: Chr 1 (100 Mb, Additive)")
simulation_h10 = load_simulation_data("./h10_simulated_ROH_Summary.parquet_dataset", "Simulated: Chr 1 (100 Mb, Dominant)")

empirical_all_autosomes_ccdf_cM = extract_official_ccdf_cM(empirical_all_autosomes)
empirical_whole_chr1_ccdf_cM = extract_official_ccdf_cM(empirical_whole_chr1)
empirical_100Mb_chr1_ccdf_cM = extract_official_ccdf_cM(empirical_100Mb_chr1)
simulation_h00_ccdf_cM = extract_official_ccdf_cM(simulation_h00)
simulation_h05_ccdf_cM = extract_official_ccdf_cM(simulation_h05)
simulation_h10_ccdf_cM = extract_official_ccdf_cM(simulation_h10)

all_plot_data_cM = bind_rows(
  empirical_all_autosomes_ccdf_cM,
  empirical_whole_chr1_ccdf_cM,
  empirical_100Mb_chr1_ccdf_cM,
  simulation_h00_ccdf_cM,
  simulation_h05_ccdf_cM,
  simulation_h10_ccdf_cM
)

all_plot_data_cM$condition = factor(all_plot_data_cM$condition, levels = c(
  "Empirical: All Autosomes", 
  "Empirical: Whole Chr 1", 
  "Empirical: Chr1 (100 Mb)",
  "Simulated: Chr 1 (100 Mb, Recessive)",
  "Simulated: Chr 1 (100 Mb, Additive)",
  "Simulated: Chr 1 (100 Mb, Dominant)"
))

all_plot_data_cM$pop_plotting = factor(all_plot_data_cM$pop_plotting, levels = c("AFR", "EUR", "EAS"))



p_ccdf_cM = ggplot(all_plot_data_cM, aes(x = cM_length, y = ccdf_y, color = condition, linetype = condition)) +
  geom_line(linewidth = 0.8) +  
  scale_y_log10(labels = label_log()) +
  scale_x_sqrt(breaks = c(1, 5, 10, 20, 50, 100, 150))+
  facet_wrap(~ pop_plotting, ncol = 1) +
  
  scale_color_manual(values = c(
    "Empirical: All Autosomes" = "grey75",
    "Empirical: Whole Chr 1"   = "grey45",
    "Empirical: Chr1 (100 Mb)"   = "black",
    "Simulated: Chr 1 (100 Mb, Recessive)" = "#e69f00",
    "Simulated: Chr 1 (100 Mb, Additive)"  = "#009373",
    "Simulated: Chr 1 (100 Mb, Dominant)"  = "#0072b2"  
  )) +
  scale_linetype_manual(values = c(
    "Empirical: All Autosomes" = "dotted",
    "Empirical: Whole Chr 1"   = "dashed",
    "Empirical: Chr1 (100 Mb)"   = "solid",
    "Simulated: Chr 1 (100 Mb, Recessive)" = "solid",
    "Simulated: Chr 1 (100 Mb, Additive)"  = "solid",
    "Simulated: Chr 1 (100 Mb, Dominant)"  = "solid"
  )) +
  guides(
    color = guide_legend(nrow = 3, byrow = TRUE),
    linetype = guide_legend(nrow = 3, byrow = TRUE)
  ) +
  theme_bw(base_size = 12, base_family = "serif") +
  labs(title = "Complementary Cumulative Distribution of ROH Lengths (cM)",
       x = "ROH Length (cM) [Square Root Scale]",
       y = "Proportion of ROH > X (Log10 Scale)",) +
  theme(
    strip.background = element_rect(fill = "grey90"), 
    strip.text = element_text(face = "bold", size = 11),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    axis.title = element_text(face = "bold"),
    legend.key.width = unit(1.5, "cm"),
    legend.text = element_text(size = 10)
  )

ggsave(filename = "ROH_CCDF_LongTail_Comparison_cM.pdf", plot = p_ccdf_cM, width = 10, height = 8, bg = "white")



empirical_all_autosomes_ccdf_Mb = extract_official_ccdf_Mb(empirical_all_autosomes)
empirical_whole_chr1_ccdf_Mb = extract_official_ccdf_Mb(empirical_whole_chr1)
empirical_100Mb_chr1_ccdf_Mb = extract_official_ccdf_Mb(empirical_100Mb_chr1)
simulation_h00_ccdf_Mb = extract_official_ccdf_Mb(simulation_h00)
simulation_h05_ccdf_Mb = extract_official_ccdf_Mb(simulation_h05)
simulation_h10_ccdf_Mb = extract_official_ccdf_Mb(simulation_h10)

all_plot_data_Mb = bind_rows(
  empirical_all_autosomes_ccdf_Mb,
  empirical_whole_chr1_ccdf_Mb,
  empirical_100Mb_chr1_ccdf_Mb,
  simulation_h00_ccdf_Mb,
  simulation_h05_ccdf_Mb,
  simulation_h10_ccdf_Mb
)

all_plot_data_Mb$condition = factor(all_plot_data_Mb$condition, levels = c(
  "Empirical: All Autosomes", 
  "Empirical: Whole Chr 1", 
  "Empirical: Chr1 (100 Mb)",
  "Simulated: Chr 1 (100 Mb, Recessive)",
  "Simulated: Chr 1 (100 Mb, Additive)",
  "Simulated: Chr 1 (100 Mb, Dominant)"
))

all_plot_data_Mb$pop_plotting = factor(all_plot_data_Mb$pop_plotting, levels = c("AFR", "EUR", "EAS"))



p_ccdf_Mb = ggplot(all_plot_data_Mb, aes(x = Mb_length, y = ccdf_y, color = condition, linetype = condition)) +
  geom_line(linewidth = 0.8) +  
  scale_y_log10(labels = label_log()) +
  scale_x_sqrt(breaks = c(1, 5, 10, 20, 50, 100, 150))+
  facet_wrap(~ pop_plotting, ncol = 1) +
  
  scale_color_manual(values = c(
    "Empirical: All Autosomes" = "grey75",
    "Empirical: Whole Chr 1"   = "grey45",
    "Empirical: Chr1 (100 Mb)"   = "black",
    "Simulated: Chr 1 (100 Mb, Recessive)" = "#e69f00",
    "Simulated: Chr 1 (100 Mb, Additive)"  = "#009373",
    "Simulated: Chr 1 (100 Mb, Dominant)"  = "#0072b2" 
  )) +
  scale_linetype_manual(values = c(
    "Empirical: All Autosomes" = "dotted",
    "Empirical: Whole Chr 1"   = "dashed",
    "Empirical: Chr1 (100 Mb)"   = "solid",
    "Simulated: Chr 1 (100 Mb, Recessive)" = "solid",
    "Simulated: Chr 1 (100 Mb, Additive)"  = "solid",
    "Simulated: Chr 1 (100 Mb, Dominant)"  = "solid"
  )) +
  guides(
    color = guide_legend(nrow = 3, byrow = TRUE),
    linetype = guide_legend(nrow = 3, byrow = TRUE)
  ) +
  theme_bw(base_size = 12, base_family = "serif") +
  labs(title = "Complementary Cumulative Distribution of ROH Lengths (Mb)",
       x = "ROH Length (Mb) [Square Root Scale]",
       y = "Proportion of ROH > X (Log10 Scale)",) +
  theme(
    strip.background = element_rect(fill = "grey90"), 
    strip.text = element_text(face = "bold", size = 11),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "bottom",      
    legend.direction = "horizontal", 
    axis.title = element_text(face = "bold"),
    legend.key.width = unit(1.5, "cm"),
    legend.text = element_text(size = 10)
  )

ggsave(filename = "ROH_CCDF_LongTail_Comparison_Mb.pdf", plot = p_ccdf_Mb, width = 10, height = 8, bg = "white")
