library(tidyverse)
library(mltools)

setwd("~/MyDrive/Projects/hemisphere_fingerprinting/code/modeling/dim_reduce/")

hemi <- read_csv("data-all_outcome-hemi_results.csv") %>%
  select(-`...1`)

hemi_summary <- hemi %>%
  group_by(k, group) %>%
  nest() %>%
  mutate(
    confusion = map(data, ~as.matrix.data.frame(table(.x$hemi, .x$predicted))),
    mcc = map_dbl(data, ~mcc(pred = .x$predicted, actuals = .x$hemi))
  )

ggplot(hemi_summary, aes(x = k, y = mcc, color = group)) +
  geom_point() +
  geom_line() +
  theme_bw()

class <- read_csv("data-all_outcome-class_results.csv") %>%
  select(-`...1`)

confusion_matrix <- function(true, pred) {

  all_true <- sort(unique(true))

  results <- tibble(true = true, pred = pred) %>%
    mutate(
      true = as.factor(true),
      pred = factor(pred, levels = all_true)
    )

  cm <- results %>%
    group_by(true, pred, .drop = FALSE) %>%
    summarize(
      n = n(),
    ) %>%
    pivot_wider(names_from = pred, values_from = n) %>%
    select(all_of(all_true)) %>%
    arrange(true) %>%
    column_to_rownames("true")

  return(cm)

}

class_summary <- class %>%
  group_by(k, group) %>%
  mutate(
    true_hemi = str_extract(class, ".H$"),
    pred_hemi = str_extract(predicted, ".H$"),
    true_hand = str_remove(class, ".H$"),
    pred_hand = str_remove(predicted, ".H$")
  ) %>%
  nest() %>%
  mutate(
    mcc = map_dbl(data, ~mcc(preds = .x$predicted, actuals = .x$class)),
    mcc_hemi = map_dbl(data, ~mcc(preds = .x$pred_hemi, actuals = .x$true_hemi)),
    mcc_hand = map_dbl(data, ~mcc(preds = .x$pred_hand, actuals = .x$true_hand))
  )

ggplot(class_summary, aes(x = k, y = mcc, color = group)) +
  geom_point() +
  geom_line() +
  theme_bw()

ggplot(class_summary, aes(x = k, y = mcc_hemi, color = group)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  labs(title = "Hemisphere prediction from k features")

ggplot(class_summary, aes(x = k, y = mcc_hand, color = group)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  labs(title = "Handedness prediction from k features")
