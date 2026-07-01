library(tidyverse)
library(igraph)

# Setup ====

select <- dplyr::select

setwd("~/MyDrive/Projects/hemisphere_fingerprinting/code/")
source("enrichment/00-enrichment_functions.R")

# Load data ====

hc <- qs::qread("modeling/inputs/hemiconnectome.rds")

hc2 <- hc %>%
  select(sub, hemi, contains("_"), -age_group) %>%
  pivot_longer(-c(sub, hemi)) %>%
  group_by(hemi, name) %>%
  summarize(
    mean = mean(value)
  ) %>%
  separate_wider_delim(name, delim = "_", names = c("x", "y")) %>%
  group_by(hemi) %>%
  nest()

hc3 <- hc2 %>%
  mutate(
    data2 = map(data, ~graph_to_tibble(.x, "mean"))
  )

hc4 <- hc3 %>%
  select(hemi, data2) %>%
  unnest(data2)

plot <- ggplot(hc4, aes(x = ROI1, y = ROI2, fill = value)) +
  geom_tile() +
  scale_x_discrete(labels = NULL) +
  scale_y_discrete(labels = NULL) +
  scale_fill_gradient2(transform = "reverse") +
  facet_grid(rows = vars(hemi)) +
  coord_equal() +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(legend.position = "none", axis.ticks = element_blank())

ggsave(filename = "hemiconnectomes_plot.png", plot = plot, width = 1.5,
       height = 3, units = "in")
