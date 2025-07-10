library(tidyverse)
library(mltools)

setwd("/Users/tkmd/Library/CloudStorage/GoogleDrive-td758@georgetown.edu/My Drive/Projects/hemisphere_fingerprinting/code/modeling/")

hemi_files <- list.files("results/pctile/", pattern = ".*_outcome-hemi_.*.csv",
                           full.names = TRUE)

class_files <- list.files("results/pctile/", pattern = ".*_outcome-class_.*.csv",
                         full.names = TRUE)

columns <- read_rds("../columns_at_each_pctile.rds")
csizes <- sapply(columns, length)

size_lut <- tibble(pctile = as.numeric(names(csizes)), size = csizes)

results_hemi <- tibble(f = hemi_files) %>%
  mutate(

    data = map(f,
               ~read_csv(.x, show_col_types = FALSE,
                         col_names = c("index", "sub", "ground", "predicted"),
                         skip = 1))

  ) %>%
  unnest(data) %>%
  mutate(
    across(c(ground, predicted), ~factor(.x, levels = c("LH", "RH")))
  ) %>%
  group_by(f) %>%
  nest() %>%
  mutate(

    contingency = map(data,
                      ~as.data.frame.matrix(table(.x$ground, .x$predicted))),


    acc = map_dbl(data, ~sum(.x$ground == .x$predicted) / nrow(.x)),
    mcc = map_dbl(contingency, ~mcc(confusionM = .x)),

    f = basename(f) %>%
      str_remove(".csv")

  ) %>%
  separate_wider_delim(f, delim = "_",
                       names = c("method", "outcome", "group", "pctile")) %>%
  mutate(
    across(c(method, outcome, group), ~str_remove(.x, "^.*-")),
    pctile = as.numeric(str_remove(pctile, "data-p"))
  ) %>%
  left_join(size_lut)

results_class <- tibble(f = class_files) %>%
  mutate(

    data = map(f,
               ~read_csv(.x, show_col_types = FALSE,
                         col_names = c("index", "sub", "ground", "predicted"),
                         skip = 1))

  ) %>%
  unnest(data) %>%
  mutate(
    across(c(ground, predicted),
           ~factor(.x, levels = c("leftyLH", "leftyRH", "rightyLH", "rightyRH")))
  ) %>%
  group_by(f) %>%
  nest() %>%
  mutate(

    contingency = map(data,
                      ~as.data.frame.matrix(table(.x$ground, .x$predicted))),


    acc = map_dbl(data, ~sum(.x$ground == .x$predicted) / nrow(.x)),
    mcc = map_dbl(contingency, ~mcc(confusionM = .x)),

    f = basename(f) %>%
      str_remove(".csv")

  ) %>%
  separate_wider_delim(f, delim = "_",
                       names = c("method", "outcome", "group", "pctile")) %>%
  mutate(
    across(c(method, outcome, group), ~str_remove(.x, "^.*-")),
    pctile = as.numeric(str_remove(pctile, "data-p"))
  ) %>%
  left_join(size_lut)

ggplot(filter(results_hemi, pctile > 0.6), aes(x = size, y = mcc)) +
  geom_point(aes(color = group)) +
  geom_line(aes(color = group)) +
  scale_x_log10(limits = c(5, NA),
                breaks = c(5, 100, 1000, 10000),
                sec.axis = sec_axis(~ .x / 16100 * 100,
                                    name = "% of all edges",
                                    breaks = c(0.1, 1, 10, 50),
                                    labels = c("0.1%", "1%", "10%", "50%"))) +
  facet_wrap(vars(method)) +
  theme_bw() +
  labs(x = "# edges", y = "MCC", title = "Two-way hemisphere: MCC",
       caption = "Omitted pctile=50,60% as they were all at ceiling") +
  theme(axis.text.x.bottom = element_text(angle = 45, hjust = 1))

ggplot(filter(results_hemi, pctile > 0.6), aes(x = size, y = acc)) +
  geom_point(aes(color = group)) +
  geom_line(aes(color = group)) +
  scale_x_log10(limits = c(5, NA),
                breaks = c(5, 100, 1000, 10000),
                sec.axis = sec_axis(~ .x / 16100 * 100,
                                    name = "% of all edges",
                                    breaks = c(0.1, 1, 10, 50),
                                    labels = c("0.1%", "1%", "10%", "50%"))) +
  scale_y_continuous(limits = c(0.8, NA)) +
  facet_wrap(vars(method)) +
  theme_bw() +
  labs(x = "# edges", y = "Accuracy", title = "Two-way hemisphere: Accuracy",
       caption = "Omitted pctile=50,60% as they were all at ceiling") +
  theme(axis.text.x.bottom = element_text(angle = 45, hjust = 1))

ggplot(filter(results_class), aes(x = size, y = mcc)) +
  geom_point(aes(color = group)) +
  geom_line(aes(color = group)) +
  scale_x_log10(limits = c(5, NA),
                breaks = c(5, 25, 100, 1000, 10000),
                sec.axis = sec_axis(~ .x / 16100 * 100,
                                    name = "% of all edges",
                                    breaks = c(0.1, 1, 10, 50),
                                    labels = c("0.1%", "1%", "10%", "50%"))) +
  facet_wrap(vars(method)) +
  theme_bw() +
  labs(x = "# edges", y = "MCC", title = "Four-way class: MCC",
       caption = "Omitted pctile=50,60% as they were all at ceiling") +
  theme(axis.text.x.bottom = element_text(angle = 45, hjust = 1))

ggplot(filter(results_class), aes(x = size, y = acc)) +
  geom_point(aes(color = group)) +
  geom_line(aes(color = group)) +
  scale_x_log10(limits = c(5, NA),
                breaks = c(5, 25, 100, 1000, 10000),
                sec.axis = sec_axis(~ .x / 16100 * 100,
                                    name = "% of all edges",
                                    breaks = c(0.1, 1, 10, 50),
                                    labels = c("0.1%", "1%", "10%", "50%"))) +
  facet_wrap(vars(method)) +
  theme_bw() +
  labs(x = "# edges", y = "MCC", title = "Four-way class: Accuracy",
       caption = "Omitted pctile=50,60% as they were all at ceiling") +
  theme(axis.text.x.bottom = element_text(angle = 45, hjust = 1))

