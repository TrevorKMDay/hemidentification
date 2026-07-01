library(tidyverse)

# Setup ====

setwd(str_glue("~/Google Drive/My Drive/Projects/hemisphere_fingerprinting/",
               "code/modeling_new/06-handedness"))
source("../../bootstrap_mcc.R")

library(parallel)
library(doParallel)
library(foreach)

cores <- detectCores()
cl <- makeCluster(cores[1] - 1) #not to overload your computer
registerDoParallel(cl)

# Load data =====

files <- list.files("results/", pattern = "*.csv", full.names = TRUE) %>%
  str_subset("hco|lang", negate = TRUE)

folds <- read_csv("../folds_to_use2.csv", show_col_types = FALSE) %>%
  select(sub, group)

data <- tibble(f = files) %>%
  mutate(
    bn = basename(f),
    data = map(f, read_csv, show_col_types = FALSE)
  ) %>%
  separate_wider_delim(bn, delim = "_", 
                       names = c("input", "clf", "outcome", NA)) %>%
  select(input, clf, outcome, data) %>%
  unnest(data) %>%
  mutate(

    outcome = str_remove(outcome, "outcome-"),
    gender_comparison = str_detect(outcome, "gender"),

    ground = recode_values(
      outcome,
      "handedness" ~ handedness,
      "handedness.gender" ~ handedness.gender,
      "class" ~ class,
      "class.gender" ~ class.gender
    ),

    ground_hemi = str_extract(ground, "[LR]H"),
    ground_gender = str_extract(ground, "[.][MF]") %>%
      str_remove("[.]"),
    ground_hand = str_remove_all(ground, "[LR]H|[.][MF]"),

    pred_hemi = str_extract(predicted, "[LR]H"),
    pred_gender = str_extract(predicted, "[.][MF]") %>%
      str_remove("[.]"),
    pred_hand = str_remove_all(predicted, "[LR]H|[.][MF]"),

  ) %>%
  select(input, clf, outcome, gender_comparison, 
          sub, gender, starts_with("ground"), 
          starts_with("pred")) %>%
  left_join(folds)

mcc_by_fold <- data %>%
  group_by(input, clf, outcome, group, gender_comparison) %>%
  nest() 

hand_only <- mcc_by_fold %>%
  filter(
    outcome == "handedness"
  ) %>%
  mutate(

    # Bootstrap MCC
    mcc_overall = map(data, bootstrap_mcc,
                      .progress = TRUE),

    mcc_hand = map(data, bootstrap_mcc, predicted = "pred_hand", 
                        ground = "ground_hand",
                      .progress = TRUE),

  ) %>%
  unnest(starts_with("mcc_"), names_sep = ".")

class <- mcc_by_fold %>%
  filter(
    outcome == "class"
  ) %>%
  mutate(

    # Bootstrap MCC
    mcc_overall = map(data, bootstrap_mcc,
                      .progress = TRUE),

    mcc_hemi = map(data, bootstrap_mcc, predicted = "pred_hemi", 
                        ground = "ground_hemi",
                      .progress = TRUE),

    mcc_hand = map(data, bootstrap_mcc, predicted = "pred_hand", 
                        ground = "ground_hand",
                      .progress = TRUE),

  ) %>%
  unnest(starts_with("mcc_"), names_sep = ".")

hand_gender <- mcc_by_fold %>%
  filter(
    outcome == "handedness.gender"
  ) %>%
  mutate(

    # Bootstrap MCC
    mcc_overall = map(data, bootstrap_mcc,
                      .progress = TRUE),

    mcc_hand = map(data, bootstrap_mcc, predicted = "pred_hand", 
                        ground = "ground_hand",
                      .progress = TRUE),

    mcc_gender = map(data, bootstrap_mcc, predicted = "pred_gender", 
                        ground = "ground_gender",
                      .progress = TRUE),

  ) %>%
  unnest(starts_with("mcc_"), names_sep = ".")

class_gender <- mcc_by_fold %>%
  filter(
    outcome == "class.gender"
  ) %>%
  mutate(

    # Bootstrap MCC
    mcc_overall = map(data, bootstrap_mcc,
                      .progress = TRUE),

    mcc_hemi = map(data, bootstrap_mcc, predicted = "pred_hemi", 
                        ground = "ground_hemi",
                      .progress = TRUE),

    mcc_hand = map(data, bootstrap_mcc, predicted = "pred_hand", 
                        ground = "ground_hand",
                      .progress = TRUE),

    mcc_gender = map(data, bootstrap_mcc, predicted = "pred_gender", 
                        ground = "ground_gender",
                      .progress = TRUE),

  ) %>%
  unnest(starts_with("mcc_"), names_sep = ".")

# Summary =====

mcc_overall <- data %>%
  group_by(input, clf, outcome, gender_comparison) %>%
  nest() 

hand_only_summary <- mcc_overall %>%
  filter(
    outcome == "handedness"
  ) %>%
  mutate(

    # Bootstrap MCC
    mcc_overall = map(data, bootstrap_mcc,
                      .progress = TRUE),

    mcc_hand = map(data, bootstrap_mcc, predicted = "pred_hand", 
                        ground = "ground_hand",
                      .progress = TRUE),

  ) %>%
  unnest(starts_with("mcc_"), names_sep = ".")

class_summary <- mcc_overall %>%
  filter(
    outcome == "class"
  ) %>%
  mutate(

    # Bootstrap MCC
    mcc_overall = map(data, bootstrap_mcc,
                      .progress = TRUE),

    mcc_hemi = map(data, bootstrap_mcc, predicted = "pred_hemi", 
                        ground = "ground_hemi",
                      .progress = TRUE),

    mcc_hand = map(data, bootstrap_mcc, predicted = "pred_hand", 
                        ground = "ground_hand",
                      .progress = TRUE),

  ) %>%
  unnest(starts_with("mcc_"), names_sep = ".")

hand_gender_summary <- mcc_overall %>%
  filter(
    outcome == "handedness.gender"
  ) %>%
  mutate(

    # Bootstrap MCC
    mcc_overall = map(data, bootstrap_mcc,
                      .progress = TRUE),

    mcc_hand = map(data, bootstrap_mcc, predicted = "pred_hand", 
                        ground = "ground_hand",
                      .progress = TRUE),

    mcc_gender = map(data, bootstrap_mcc, predicted = "pred_gender", 
                        ground = "ground_gender",
                      .progress = TRUE),

  ) %>%
  unnest(starts_with("mcc_"), names_sep = ".")

class_gender_summary <- mcc_overall %>%
  filter(
    outcome == "class.gender"
  ) %>%
  mutate(

    # Bootstrap MCC
    mcc_overall = map(data, bootstrap_mcc,
                      .progress = TRUE),

    mcc_hemi = map(data, bootstrap_mcc, predicted = "pred_hemi", 
                        ground = "ground_hemi",
                      .progress = TRUE),

    mcc_hand = map(data, bootstrap_mcc, predicted = "pred_hand", 
                        ground = "ground_hand",
                      .progress = TRUE),

    mcc_gender = map(data, bootstrap_mcc, predicted = "pred_gender", 
                        ground = "ground_gender",
                      .progress = TRUE),

  ) %>%
  unnest(starts_with("mcc_"), names_sep = ".")

# Plotting =====

results <- bind_rows(hand_only, class, hand_gender, class_gender) %>%
  pivot_longer(starts_with("mcc_")) %>%
  separate_wider_delim(name, delim = ".", names = c("mcc", "x")) %>%
  pivot_wider(names_from = x)

results_summary <- bind_rows(hand_only_summary, class_summary, 
                             hand_gender_summary, class_gender_summary) %>%
  pivot_longer(starts_with("mcc_")) %>%
  separate_wider_delim(name, delim = ".", names = c("mcc", "x")) %>%
  pivot_wider(names_from = x)

results_lda <- results %>%
  filter(
    clf == "clf-lda"
  ) %>%
  mutate(
    mcc = str_remove(mcc, "mcc_") %>%
      str_to_sentence(),
    input = toupper(str_remove(input, "data-"))
  )

results_summary_lda <- results_summary %>%
  filter(
    clf == "clf-lda"
  ) %>%
  mutate(
    mcc = str_remove(mcc, "mcc_") %>%
      str_to_sentence(),
    input = toupper(str_remove(input, "data-"))
  )

ggplot(results_lda, aes(x = gender_comparison, color = mcc)) +
  geom_pointrange(
    aes(y = mcc_mean, ymin = mcc_ci_lo, ymax = mcc_ci_hi),
    position = position_dodge2(width = .5)
  ) +
  geom_line(
    data = filter(results_summary_lda, mcc == "Hand"), 
    aes(x = gender_comparison, y = mcc_mean, group = input),
    linewidth = 1.5
  ) +
  geom_pointrange(
    data = results_summary_lda, 
    aes(y = mcc_mean, ymin = mcc_ci_lo, ymax = mcc_ci_hi, fill = mcc),
    position = position_dodge2(width = .5),
    shape = 21, color = "black", size = 0.9
  ) +
  geom_hline(yintercept = 0, color = "red") +
  scale_x_discrete(labels = c("No", "Yes")) +
  facet_wrap(vars(input), scales = "free_x") +
  theme_bw() + 
  labs(x = "With gender", y = "MCC (95% CI)", color = "Prediction", 
        fill = "Prediction") +
  theme(legend.position = "bottom")

ggsave("plots/gender_changeinMCC.png", width = 6.5, height = 4)

gender_effect <- aov(mcc_mean ~ group + input + gender_comparison,
                      data = filter(results_lda, mcc == "mcc_hand"))

summary(gender_effect)

range(results_lda$mcc_mean[results_lda$mcc == "mcc_gender"], na.rm = TRUE)
