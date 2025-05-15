library(tidyverse)
library(mltools)
library(doParallel)
library(foreach)

# Parallel setup

cores <- detectCores()
cl <- makeCluster(cores[1] - 1) #not to overload your computer
registerDoParallel(cl)

# Function to boostrap mcc

bootstrap_mcc <- function(data, R = 1000, ground = "ground",
                          predicted = "predicted", levels = NULL) {

  estimates <- rep(NA, R)

  estimates <- foreach (i = 1:R, .combine = "c") %dopar% {

    new_data <- dplyr::slice_sample(data, prop = 1, replace = TRUE)

    g <- unlist(new_data[, ground], use.names = FALSE)
    p <- unlist(new_data[, predicted], use.names = FALSE)

    if (is.null(levels))
      levels <- unique(g)

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


# bootstrap_mcc(results2$data[[1]])

setwd("~/Google Drive/My Drive/Projects/hemisphere_fingerprinting/code/modeling")

files_base <- list.files("results/base/", pattern = "*.csv", full.names = TRUE)

files_trans <- list.files("results/trans/", pattern = "*.csv",
                          full.names = TRUE)

files_osampl <- list.files("results/osampl/", pattern = "*.csv",
                          full.names = TRUE)

files_langtask <- list.files("results/langtask/", pattern = "*.csv",
                             full.names = TRUE)

files_full <- list.files("results/full/", pattern = "*.csv", full.names = TRUE)

all_files1 <- c(files_base, files_trans, files_osampl, files_langtask,
                files_full)

# files_nets <- list.files("results/networks/", pattern = "*.csv",
                         # full.names = TRUE)

results <- tibble(f = all_files1) %>%
  mutate(
    data = map(f, read_csv, show_col_types = FALSE,
               name_repair = "unique_quiet", .progress = TRUE),
    f = basename(f) %>%
      str_remove(".csv")
  )

results %>%
  filter(
    str_detect(f, "data-base"),
    str_detect(f, "hands-all")
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
           ~str_remove(.x, "^.*-"))

  ) %>%
  unnest(data) %>%
  mutate(
    ground_hand = str_remove(ground, ".H$"),
    ground_hemi = str_remove(predicted, "^(left|right)y"),
    pred_hand = str_remove(predicted, ".H$"),
    pred_hemi = str_remove(predicted, "^(left|right)y"),
    across(c(starts_with("ground"), starts_with("pred")), as.factor),
  ) %>%
  select(-`...1`, -f) %>%
  group_by(outcome, method, group, hands, hemis, input, n)

results_big %>%
  filter(
    !(ground_hemi %in% c("LH", "RH"))
  )

results2 <- results_big %>%
  filter(
    input %in% c("base", "osampl", "langtask", "full")
  ) %>%
  nest() %>%
  select(outcome, method, group, hands, hemis, input, data, n)

results_trans <-  results_big %>%
  filter(
    input %in% c("trans")
  ) %>%
  nest() %>%
  select(outcome, method, group, hands, hemis, input, data, n)

## Results by outcome ====

### Predicting handedness from hemisphere connectivity ====

results_hands <- results2 %>%
  filter(
    hands == "all",
    outcome %in% c("hand", "hand2")
  ) %>%
  mutate(

    accuracy = map_dbl(data, ~sum(.x$ground_hand == .x$pred_hand) / nrow(.x)),

    # Specify only evaluating handedness
    # WARNING: DON'T KILL THIS
    mcc = map(data, bootstrap_mcc, R = 100,
              ground = "ground_hand", predicted = "pred_hand",
              .progress = TRUE)

  ) %>%
  unnest(mcc)

results_hands %>%
  filter(
    input == "base"
  ) %>%
  pull(mcc_mean) %>%
  range()

results_hands %>%
  filter(
    input == "osampl"
  ) %>%
  pull(mcc_mean) %>%
  range()

results_hands %>%
  filter(
    input == "langtask"
  ) %>%
  pull(mcc_mean) %>%
  range()

results_hands %>%
  filter(
    input == "full"
  ) %>%
  pull(mcc_mean) %>%
  range()

results_hands %>%
  filter(
    method == "lda",
    hands == "all",
    input == "base",
    hemis == "all"
  )

summary_hands <- results_hands %>%
  filter(
    outcome %in% c("hand", "hand2")
  ) %>%
  group_by(method, outcome, hands, input, hemis) %>%
  summarize(
    n2 = n(),
    mcc_hands_mean = weighted.mean(mcc_mean, n),
    mcc_hands_sd = sqrt(Hmisc::wtd.var(mcc_mean, n))
  ) %>%
  mutate(
    mcc_hands_se = mcc_hands_sd / sqrt(n2)
  )

mcc_quick_model <- results_hands %>%
  ungroup() %>%
  filter(
    hemis == "all"
  ) %>%
  select(-hemis, -data) %>%
  fastDummies::dummy_cols(c("method", "group"))

mcc_lm <- aov(mcc_mean ~ method + group,  data = mcc_quick_model)

TukeyHSD(mcc_lm)

png("plots/predict-hand_only.png", width = 5.5, height = 4, units = "in",
    res = 300)

ggplot(results_hands, aes(x = interaction(input, method))) +
  geom_pointrange(aes(y = mcc_mean, ymin = mcc_ci_lo, ymax = mcc_ci_hi,
                      group = group),
                  position = position_dodge(width = 0.5),
                  size = 0.25, alpha = 0.25) +
  geom_pointrange(data = summary_hands,
                  aes(x = interaction(input, method), y = mcc_hands_mean,
                      ymin = mcc_hands_mean - mcc_hands_se,
                      ymax = mcc_hands_mean + mcc_hands_se,
                      shape = input, linetype = input,
                      color = input),
                  size = 0.5) +
  geom_hline(yintercept = 0, color = "red") +
  scale_y_continuous(limits = c(NA, NA)) +
  facet_grid(cols = vars(hemis), rows = vars(outcome), scales = "free_x") +
  theme_bw() +
  labs(x = "Method", y = "MCC", title = "Predict handedness (only)",
       color = "Test group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

dev.off()

results_hands %>%
  filter(
    method == "lda",
    hemis == "all"
  ) %>%
  distinct()

### Predicting hemisphere ====

results_hemis <- results2 %>%
  filter(
    hemis == "all",
    input != "trans",
    outcome == "hemi"
  ) %>%
  mutate(

    # Specify only evaluating hemi
    mcc = map(data, bootstrap_mcc, R = 1000,
              ground = "ground_hemi", predicted = "pred_hemi",
              levels = c("LH", "RH"))

  ) %>%
  unnest(mcc)

summary_hemis <- results_hemis %>%
  group_by(method, hemis, hands, input) %>%
  summarize(
    n2 = n(),
    mcc_hemis_mean = weighted.mean(mcc, n),
    mcc_hemis_sd = sqrt(Hmisc::wtd.var(mcc, n))
  ) %>%
  mutate(
    mcc_hemis_se = mcc_hemis_sd / sqrt(n2)
  )

ggplot(results_hemis, aes(x = interaction(input, method))) +
  geom_pointrange(aes(y = mcc_mean, ymin = mcc_ci_lo, ymax = mcc_ci_hi,
                      color = group),
                  position = position_dodge(width = 0.5),
                  size = 0.25, alpha = 0.75) +
  geom_pointrange(data = summary_hemis,
                  aes(x = interaction(input, method), y = mcc_hemis_mean,
                      ymin = mcc_hemis_mean - mcc_hemis_se,
                      ymax = mcc_hemis_mean + mcc_hemis_se,
                      shape = input, linetype = input),
                  size = 0.5) +
  geom_hline(yintercept = 0, color = "red") +
  scale_y_continuous(limits = c(NA, NA)) +
  facet_grid(cols = vars(hands)) +
  theme_bw() +
  labs(x = "Method", y = "MCC", title = "Predict two-way hemisphere (only)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

## Predicting four-way handedness from hemisphere connectivity

zero_rows <- function(confusion, rows) {

  x <- confusion
  x[rows, ] <- 0

  return(x)

}

four_levels <- c("leftyLH", "leftyRH", "rightyLH", "rightyRH")

# This is the model trained on all at once
results_hemis4way <- results2 %>%
  filter(
    hemis == "all",
    hands == "all",
    input != "trans"
  ) %>%
  mutate(

    lefties = map(data, ~filter(.x, ground_hand == "lefty")),
    righties = map(data, ~filter(.x, ground_hand == "righty")),

    all = map(data, bootstrap_mcc, R = 1000, levels = four_levels),

    lefties = map(lefties, bootstrap_mcc, R = 1000, levels = four_levels),
    righties = map(righties, bootstrap_mcc, R = 1000, levels = four_levels),

  ) %>%
  unnest(c(all, lefties, righties), names_sep = "_")

results_hemis4way2 <- results_hemis4way %>%
  ungroup() %>%
  select(method, input, group, n, contains("mcc")) %>%
  pivot_longer(-c(method, input, group, n)) %>%
  separate_wider_delim(name, delim = "_", names = c("hgroup", "name"),
                       too_many = "merge") %>%
  pivot_wider()

summary_hemis4way <- results_hemis4way2 %>%
  group_by(method, input, hgroup) %>%
  summarize(
    n2 = n(),
    mcc_mean = weighted.mean(mcc_mean, n),
    mcc_sd = sqrt(Hmisc::wtd.var(mcc_mean, n)),
  ) %>%
  mutate(
    mcc_se = mcc_sd / sqrt(n2)
  )

png("plots/predict-4way.png", width = 5.5, height = 4, units = "in",
    res = 300)

ggplot(results_hemis4way2, aes(x = interaction(input, method))) +
  geom_pointrange(aes(y = mcc_mean, ymin = mcc_ci_lo, ymax = mcc_ci_hi,
                      color = group),
                  position = position_dodge(width = 0.5),
                  size = 0.25, alpha = 0.75) +
  geom_pointrange(data = summary_hemis4way,
                  aes(x = interaction(input, method),
                      y = mcc_mean,
                      ymin = mcc_mean - mcc_se,
                      ymax = mcc_mean + mcc_se,
                      shape = input, linetype = input),
                  size = 0.5) +
  geom_hline(yintercept = 0, color = "red") +
  coord_cartesian(ylim = c(NA, 1)) +
  facet_wrap(vars(hgroup)) +
  theme_bw() +
  labs(x = "Method", y = "MCC", title = "Predict four-way hemisphere") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

dev.off()

# Transconnectome ====

trans_big <- results_trans %>%
  unnest(data) %>%
  ungroup() %>%
  select(-contains("hemi"), -hands, -hemis, -outcome, -input, -ground,
         -predicted) %>%
  group_by(method, group) %>%
  nest() %>%
  mutate(

    confusion = map(data, ~as.data.frame.matrix(table(.x$ground_hand,
                                                      .x$pred_hand))),

    mcc = map_dbl(confusion, ~mcc(confusionM = .x))

  )

png("plots/predict-hand_trans.png", width = 5.5, height = 4, units = "in",
    res = 300)

ggplot(trans_big, aes(x = method, y = mcc, color = group)) +
  geom_point() +
  theme_bw() +
  labs(title = "Predict handedness from trans connections")

dev.off()


