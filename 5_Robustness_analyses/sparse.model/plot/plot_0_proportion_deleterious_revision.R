library(tidyverse)
library(ggplot2) 
library(dplyr)
options(scipen = 999)
env_list = c('1.0')
x_list = c('1.0', '0.3')

for (x in x_list){
  for (env in env_list){
    df_raw = read_csv(paste0('summary_different_kind_of_inds.subset.',x,'.csv'))
    df_raw = df_raw %>%
      rename(ROH_type = category) %>%
      pivot_longer(
        cols = c(`len>0_score!=0`, `len>0_score==0`, `len==0_score==0`), 
        names_to = "Status", 
        values_to = "Count"
      ) %>%
      mutate(
        hXX = case_when(
          hXX == "h00" ~ "recessive",
          hXX == "h05" ~ "additive",
          hXX == "h10" ~ "dominant",
          hXX == "h15" ~ "mixed",
          TRUE ~ hXX 
        ),
        ROH_type = case_when(
          ROH_type == "A" ~ "Short ROH",
          ROH_type == "B" ~ "Medium ROH",
          ROH_type == "C" ~ "Long ROH",
          ROH_type == "Non-ROH" ~ "Non-ROH",
          TRUE ~ ROH_type
        ),
        Status = case_when(
          Status == "len==0_score==0" ~ "No ROH",
          Status == "len>0_score!=0" ~ "Have at least one non-zero contribution ROH",
          Status == "len>0_score==0" ~ "Only have zero contribution ROH",
          TRUE ~ Status
        ),
        pop = case_when(
          pop == "p1" ~ "African",
          pop == "p2" ~ "European",
          pop == "p3" ~ "East Asian",
          TRUE ~ pop
        ))
    df_filtered = df_raw[df_raw$ROH_type != "Non-ROH", ]
    df_summary <- df_filtered %>%
      mutate(
        #ROH_Type = factor(ROH_Type, levels = c("Short ROH", "Medium ROH", "Long ROH", "Non-ROH")),
        ROH_type = factor(ROH_type, levels = c("Short ROH", "Medium ROH", "Long ROH")),
        pop = factor(pop, levels = c("African", "European", "East Asian")),
        hXX = factor(hXX, levels = c("recessive", "additive", "dominant", "mixed")),
        Status = factor(Status, levels = c("No ROH", 
                                           "Only have zero contribution ROH", 
                                           "Have at least one non-zero contribution ROH"))
      )
    
    # Calculate the proportion instead of just counting numbers
    df_summary$Ratio = df_summary$Count/500
    df_summary = df_summary %>% select(-"Count")
    
    # Calculate the mean value and standard error of each group
    df_ready = df_summary %>% group_by(hXX,pop,ROH_type,Status) %>%
      summarise(mean_Ratio = mean(Ratio),
                sd_Ratio = sd(Ratio),
                .groups = 'drop')
    df_ready$se_Ratio = df_ready$sd_Ratio/sqrt(500)
    
    
    df_plot = ggplot(df_ready, aes(x = ROH_type, y = mean_Ratio, fill = Status)) +
      geom_col(position = "stack", alpha = 0.9, width = 0.65) +
      scale_fill_manual(values = c(
        "No ROH" = rgb(148, 203, 236, maxColorValue = 255),
        "Only have zero contribution ROH" = rgb(220, 205, 125, maxColorValue = 255),
        "Have at least one non-zero contribution ROH" = rgb(194, 106, 119, maxColorValue = 255)
      ))+
      facet_grid(hXX ~ pop) +
      labs(
        title = bquote(bold("x = " * .(x))),
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
        legend.direction = "horizontal",
        axis.title = element_text(face = "bold"), 
        legend.title = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5)
      )
    df_plot
    ggsave(
      filename = paste0('zero_ind_proportion_deleterious.subset.',x,'.pdf'),
      plot     = df_plot,
      device   = cairo_pdf, 
      width    = 8,
      height   = 8,
      units    = "in"
    )
  }
}