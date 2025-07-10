library(tidyverse)
library(mltools)
library(doParallel)
library(foreach)

select <- dplyr::select

# Parallel setup ====

cores <- detectCores()
cl <- makeCluster(cores[1] - 1) #not to overload your computer
registerDoParallel(cl)

# Function to boostrap MCC ====

inputs0 <- tibble(f = list.files("inputs/subsets/", full.names = TRUE)) %>%
  mutate(
    data = map(f, reticulate::py_load_object),
  )

inputs <- inputs0 %>%
  mutate(
    ncol = map_int(data, ncol),
    pctile = str_remove_all(f, "inputs/subsets//hemiconnectome_pctile|.pickle") %>%
      as.numeric()
  ) %>%
  select(-f, -data) %>%
  arrange(pctile)

bootstrap_mcc <- function(data, R = 1000, ground = "ground",
                          predicted = "predicted", levels = NULL) {

  estimates <- rep(NA, R)

  estimates <- foreach (i = 1:R, .combine = "c") %dopar% {

    new_data <- dplyr::slice_sample(data, prop = 1, replace = TRUE)

    g <- unlist(new_data[, ground], use.names = FALSE)
    p <- unlist(new_data[, predicted], use.names = FALSE)

    values <- unique(c(g, p))

    if (!is.null(levels))
      if (length(values) != length(levels))
        warning("Levels do not match values in data")
    else
      levels <- values

    g <- factor(g, levels = levels)
    p <- factor(p, levels = levels)

    conf <- as.data.frame.matrix(table(g, p))
    conf2 <- conf[levels, levels]

    mltools::mcc(confusionM = conf2)

  }

  est_mean <- mean(estimates)
  est_sd <- sd(estimates)
  est_ci <- est_sd * 2

  return(tibble(mcc_mean = est_mean,
                mcc_ci_lo = est_mean - est_ci,
                mcc_ci_hi = est_mean + est_ci))

}

# Setup ====

setwd("~/Google Drive/My Drive/Projects/hemisphere_fingerprinting/code/modeling")

# Find data ====

files_pctile <- list.files("results/pctile//", pattern = "*.csv",
                           full.names = TRUE)


results <- tibble(f = files_pctile) %>%
  mutate(
    data = map(f, read_csv, show_col_types = FALSE,
               name_repair = "unique_quiet", .progress = TRUE),
    f = basename(f) %>%
      str_remove(".csv")
  )

results_big <- results %>%
  mutate(

    n = map_int(data, nrow),

    method = str_extract(f, "method-[^_]*"),
    outcome = str_extract(f, "outcome-[^_]*"),
    group = str_extract(f, "test-[^_]*"),
    input = str_extract(f, "data-.*"),


    hands = replace_na(str_extract(f, "hands-[^_]*"), "all"),
    hemis = replace_na(str_extract(f, "hemi-[^_]*"), "all"),

    across(c(outcome, method, group, hands, hemis, input),
           ~str_remove(.x, "^.*-")),

    pctile = as.numeric(str_remove(input, "pctile")),

  ) %>%
  unnest(data) %>%
  mutate(
    ground_hand = str_remove(ground, ".H$"),
    ground_hemi = str_remove(ground, "^(left|right)y"),
    pred_hand = str_remove(predicted, ".H$"),
    pred_hemi = str_remove(predicted, "^(left|right)y"),
  ) %>%
  select(-`...1`, -f) %>%
  group_by(outcome, method, group, hands, hemis, input, pctile, n)

results2 <- results_big %>%
  nest() %>%
  mutate(

    # Specify only evaluating hemi
    all = map(data, bootstrap_mcc, R = 100,
              ground = "ground", predicted = "predicted",
              levels = c("leftyLH", "leftyRH", "rightyLH", "rightyRH")),

    hemi = map(data, bootstrap_mcc, R = 100,
                   ground = "ground_hemi", predicted = "pred_hemi",
                   levels = c("LH", "RH")),

    hand = map(data, bootstrap_mcc, R = 100,
               ground = "ground_hand", predicted = "pred_hand",
               levels = c("lefty", "righty")),

  ) %>%
  unnest(c(all, hemi, hand), names_sep = "_") %>%
  left_join(inputs)

my_breaks <- c(25, 200, 800, 4000)

ggplot(results2, aes(x = ncol, y = all_mcc_mean, color = group)) +
  geom_point() +
  geom_line() +
  scale_x_log10(breaks = my_breaks) +
  scale_y_continuous(limits = c(NA, NA)) +
  facet_wrap(vars(method)) +
  theme_bw() +
  labs(title = "Four-way classification")

ggplot(results2, aes(x = ncol, y = hand_mcc_mean, color = group)) +
  geom_point() +
  geom_line() +
  scale_x_log10(breaks = my_breaks) +
  scale_y_continuous(limits = c(NA, 1)) +
  facet_wrap(vars(method)) +
  theme_bw() +
  labs(title = "Handedness classification")

ggplot(results2, aes(x = ncol, y = hemi_mcc_mean, color = group)) +
  geom_point() +
  geom_line() +
  scale_x_log10(breaks = my_breaks) +
  scale_y_continuous(limits = c(NA, NA)) +
  facet_wrap(vars(method)) +
  theme_bw() +
  labs(title = "Hemispheric classification")

# double-check my work

