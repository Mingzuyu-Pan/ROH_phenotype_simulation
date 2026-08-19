library(tidyverse)

df <- read.csv("Mb_length_in_cM_class.csv") %>%
  mutate(
    class = case_when(
      class == "A" ~ "Short ROH",
      class == "B" ~ "Medium ROH",
      class == "C" ~ "Long ROH",
      TRUE ~ class
    ),
    class = factor(class, levels = c("Short ROH", "Medium ROH", "Long ROH")),
    pop = case_when(
      pop == "p1" ~ "African",
      pop == "p2" ~ "European",
      pop == "p3" ~ "East Asian",
      TRUE ~ pop
    ),
    pop = factor(pop, levels = c("African", "European", "East Asian"))
  )
pop_colors <- c("African" = "#F4B5C1", "European" = "#F5E6A3", "East Asian" = "#A8D4E6")


p1 = ggplot(df, aes(x = class, y = length_Mb, fill = pop)) +
  geom_boxplot(linewidth = 0.3) +
  facet_wrap(~ pop) +
  scale_fill_manual(values = pop_colors)+
  labs(x = "ROH class (cM)", y = "Mean ROH length (Mb)") +
  scale_y_continuous(limits = c(0, 3.0), breaks = seq(0, 3.0, 0.5))+
  theme_bw()

p2 = ggplot(df, aes(x = class, y = length_cM, fill = pop)) +
  geom_boxplot(linewidth = 0.3) +
  facet_wrap(~ pop) +
  scale_fill_manual(values = pop_colors)+
  labs(x = "ROH class (cM)", y = "Mean ROH length (cM)") +
  scale_y_continuous(limits = c(0, 3.0), breaks = seq(0, 3.0, 0.5))+
  theme_bw()

ggsave("Mb_length_in_cM_class_boxplot.pdf", plot = p1, width = 10, height = 4)
ggsave("cM_length_in_cM_class_boxplot.pdf", plot = p2, width = 10, height = 4)
