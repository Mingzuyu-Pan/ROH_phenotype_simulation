# author: Mingzuyu Pan(mzp5919@psu.edu)

# Load libraries for following plotting and analysis
library(tidyverse)
library(ggplot2) 
library(dplyr)
options(scipen = 999)

# Read raw data
df_raw = read_csv("calculate_0_proportion_deleterious.csv")

# Double check the existence of nan value before moving forward
print(all(df_raw$A_nan == 0))
print(all(df_raw$B_nan == 0))
print(all(df_raw$C_nan == 0))
print(all(df_raw$NONE_nan == 0))

# Delete the columns ending with nan
df_ready = df_raw %>%
  select(-ends_with("_nan"))

# Change the format of dataframe
df_analysis = df_ready %>%
pivot_longer(
  cols = A_0:NONE_2,
  names_to = c("ROH_Type", "Status"),
  names_sep = "_",
  values_to = "Count" 
) 

# Calculate the proportion instead of just counting numbers
df_analysis$Ratio = df_analysis$Count/500
df_analysis = df_analysis %>% select(-"Count")

# Calculate the mean value and standard error of each group
df_summary = df_analysis %>% group_by(hXX,pop,ROH_Type,Status) %>%
                        summarise(mean_Ratio = mean(Ratio),
                                  sd_Ratio = sd(Ratio),
                                  .groups = 'drop')
df_summary$se_Ratio = df_summary$sd_Ratio/sqrt(500)

df_summary_keep = df_summary
df_summary = df_summary[df_summary$ROH_Type != "NONE", ]

df_summary = df_summary %>%
  mutate(
    hXX = case_when(
      hXX == "h00" ~ "recessive",
      hXX == "h05" ~ "additive",
      hXX == "h10" ~ "dominant",
      hXX == "h15" ~ "mixed",
      TRUE ~ hXX 
    ),
    ROH_Type = case_when(
      ROH_Type == "A" ~ "Short ROH",
      ROH_Type == "B" ~ "Medium ROH",
      ROH_Type == "C" ~ "Long ROH",
      #ROH_Type == "NONE" ~ "Non-ROH",
      TRUE ~ ROH_Type
    ),
    Status = case_when(
      Status == "0" ~ "No ROH",
      Status == "2" ~ "Have at least one non-zero contribution ROH",
      Status == "1" ~ "Only have zero contribution ROH",
      TRUE ~ Status
    ),
    pop = case_when(
      pop == "p1" ~ "African",
      pop == "p2" ~ "European",
      pop == "p3" ~ "East Asian",
      TRUE ~ pop
    ))

df_summary <- df_summary %>%
  mutate(
    #ROH_Type = factor(ROH_Type, levels = c("Short ROH", "Medium ROH", "Long ROH", "Non-ROH")),
    ROH_Type = factor(ROH_Type, levels = c("Short ROH", "Medium ROH", "Long ROH")),
    pop = factor(pop, levels = c("African", "European", "East Asian")),
    hXX = factor(hXX, levels = c("recessive", "additive", "dominant","mixed")),
    Status = factor(Status, levels = c("No ROH", 
                                       "Only have zero contribution ROH", 
                                       "Have at least one non-zero contribution ROH"))
  )


df_plot = ggplot(df_summary, aes(x = ROH_Type, y = mean_Ratio, fill = Status)) +
  geom_col(position = "stack", alpha = 0.9, width = 0.65) +
  scale_fill_manual(values = c(
    "No ROH" = rgb(148, 203, 236, maxColorValue = 255),
    "Only have zero contribution ROH" = rgb(220, 205, 125, maxColorValue = 255),
    "Have at least one non-zero contribution ROH" = rgb(194, 106, 119, maxColorValue = 255)
  ))+
  facet_grid(hXX ~ pop) +
  labs(
    x = "Class of ROH",
    y = "Ratio",
    fill = "Status"
  ) +
  theme_minimal() +
  theme(
    text = element_text(family = "Arial", size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "bottom",  # Position the legend at the bottom
    legend.direction = "horizontal"
  )
df_plot
ggsave(
  filename = "zero_ind_proportion_deleterious.pdf",
  plot     = df_plot,
  device   = cairo_pdf, 
  width    = 8,
  height   = 8,
  units    = "in"
)

