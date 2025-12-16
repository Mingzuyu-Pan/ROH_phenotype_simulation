library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(showtext)
library(svglite)
library(patchwork)


hXX_list = c('h00','h05','h10','h15')

all_summaries_list = list()

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
  path = paste0('proportion_',hXX,'_plot_neutral.csv')
  summary = read.csv(path)
  summary$hXX = hXX
  summary_plot =  summary %>% select(-c(Rho, Short.ROH.count,Short.ROH.standard.variation,Medium.ROH.count,Medium.ROH.standard.variation,Long.ROH.count,Long.ROH.standard.variation,Non.ROH.count,Non.ROH.standard.variation))
  all_summaries_list[[hXX]] = summary_plot
}

final_combined_df = dplyr::bind_rows(all_summaries_list)

write.xlsx(
  x = final_combined_df,
  file = "./Table6_proportion_explained_se_neutral.xlsx"
)

final_combined_df_long = final_combined_df %>%
  pivot_longer(
    cols = matches("\\.(mean|standard\\.error)$"), 
    names_to = c("ROH_Category", ".value"),
    names_pattern = "(.+)\\.(mean|standard\\.error)$"
  )

desired_order = c("Short.ROH", "Medium.ROH", "Long.ROH", "Non.ROH")
final_plot_data = final_combined_df_long %>%
  mutate(ROH_Category = factor(ROH_Category, levels = desired_order))

color_palette = c("h00" ='#B5134E', "h05" = '#1E88E5', "h10" = '#FFC107', "h15" = '#5AA99B')
line_styles = c("African" = "solid", "European" = "dashed", "East Asian" = "dotted")

new_facet_titles = c(
  "Short.ROH" = "Short ROH",
  "Medium.ROH" = "Medium ROH",
  "Long.ROH" = "Long ROH",
  "Non.ROH" = "Non-ROH"
)


final_plot_simplified = ggplot(final_plot_data, aes(x = Tau, y = mean, color = hXX, linetype = Pop, group = interaction(hXX, Pop))) +
  geom_errorbar(aes(ymin = mean - standard.error, ymax = mean + standard.error, color = hXX), width = 0.04, alpha = 0.8,linetype = "solid") +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1) +
  
  facet_wrap(
    ~ ROH_Category, 
    ncol = 2, 
    scales = "free_y",
    labeller = as_labeller(new_facet_titles)
  ) +
  
  scale_x_continuous(breaks = unique(final_plot_data$Tau)) +
  scale_color_manual(values = color_palette,
                     breaks = c("h00", "h05", "h10", "h15"),
                     labels = c("Recessive", "Additive", "Dominant", "Mixed") ) +
  scale_linetype_manual(values = line_styles) +
  
  labs(
    x = "Tau",
    y = "Proportion of phenotype explained",
    color = "Genetic Model",
    linetype = "Population"
  ) +
  
  theme_bw() +
  theme(
    legend.box = "vertical",
    axis.title = element_text(face = "bold"),
    text = element_text(family = "Arial"),
    legend.position = "bottom",
    strip.background = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    #strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold")
  )

print(final_plot_simplified)

ggsave(
  filename = "proportion_with_se_neutral.pdf",
  plot     = final_plot_simplified,
  device   = cairo_pdf,  
  width    = 6,
  height   = 9,
  units    = "in"
)
