library(tidyverse)
library(mltools)
library(doParallel)
library(foreach)

select <- dplyr::select

# Parallel setup ====

cores <- detectCores()
cl <- makeCluster(cores[1] - 1) # not to overload your computer
registerDoParallel(cl)

# Setup ====

setwd(paste0("~/Google Drive/My Drive/Projects/hemisphere_fingerprinting/code/",
             "modeling/analyze_results/"))

source("bootstrap_mcc.R")

# Find data ====

all_files1 <- list.files("../results/baseZ/", pattern = ".*_data-.*.csv",
                         full.names = TRUE, recursive = TRUE)

# Read data in from files
results <- tibble(f = all_files1) %>%
  mutate(
    data = map(f, read_csv, show_col_types = FALSE,
               name_repair = "unique_quiet",
               .progress = TRUE),
    f = basename(f) %>%
      str_remove(".csv")
  )

results_big <- results %>%
  mutate(

    n = map_int(data, nrow),

    method = str_extract(f, "method-[^_]*"),
    group = str_extract(f, "test-[^_]*"),
    input = str_extract(f, "data-[^_]*"),

    subset = str_extract(f, "subset-[^_]*"),

    across(c(method, group, input, subset), ~str_remove(.x, "^.*-"))

  ) %>%
  unnest(data) %>%
  mutate(
    ground_hemi = str_extract(ground, ".H$"),
    pred_hemi = str_extract(predicted, ".H$"),
  ) %>%
  group_by(method, group, input, subset, n) %>%
  filter(
    ground %in% c("rightyRH", "rightyLH", "RH", "LH")
  ) %>%
  select(-`...1`, -f)

results_nested <- results %>%
  mutate(

    n = map_int(data, nrow),

    method = str_extract(f, "method-[^_]*"),
    group = str_extract(f, "test-[^_]*"),
    input = str_extract(f, "data-[^_]*"),
    hand = str_extract(f, "hands-[^_]*"),

    subset = str_extract(f, "subset-[^_]*") %>%
      replace_na("none"),

    across(c(method, group, input, subset, hand), ~str_remove(.x, "^.*-"))

  ) %>%
  unnest(data) %>%
  mutate(
    ground_hemi = str_extract(ground, ".H$"),
    pred_hemi = str_extract(predicted, ".H$"),
  ) %>%
  group_by(method, group, hand, input, subset, n) %>%
  filter(
    ground %in% c("rightyRH", "rightyLH", "LH", "RH"),
    subset == "none"
  ) %>%
  select(-`...1`, -f)

results_mcc <- results_nested %>%
  nest() %>%
  mutate(

    # Specify only evaluating hemi
    mcc = map(data, bootstrap_mcc,
              R = 100, ground = "ground_hemi", predicted = "pred_hemi",
              levels = c("LH", "RH"), method = 2,
              .progress = TRUE),

    acc = map_dbl(data, ~sum(.$ground_hemi == .$pred_hemi) / nrow(.)),

    n_errors = map_int(data, ~sum(.$ground_hemi != .$pred_hemi))

  ) %>%
  unnest(mcc)

results_mcc_summary <- results_mcc %>%
  ungroup() %>%
  group_by(method, input) %>%
  summarize(
    mcc_min = min(mcc_mean),
    mcc_max = max(mcc_mean),
    acc_min = min(acc),
    acc_max = max(acc),
  ) %>%
  mutate(
    across(where(is.numeric), ~round(., 3))
  )

ggplot(results_mcc, aes(x = input, y = mcc_mean, color = group)) +
  geom_pointrange(aes(ymin = mcc_ci_lo, ymax = mcc_ci_hi),
                  position = position_dodge2(width = 0.5),
                  alpha = 0.75) +
  facet_grid(cols = vars(method), rows = vars(hand)) +
  theme_bw() +
  labs(x = "Input data", y = "MCC (95% CI)")
