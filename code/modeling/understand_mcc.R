library(tidyverse)
library(mltools)

setwd("~/MyDrive/Projects/hemisphere_fingerprinting/code")

accuracy <- function(matrix) {

  correct <- sum(diag(matrix))
  total <- sum(matrix)

  return(correct / total)

}

testing_mcc <- tibble(n_lefty_correct = 0:20) %>%
  mutate(
    mat = map(n_lefty_correct,
              ~matrix(data = c(.x, 0, 20 - .x, 200), nrow = 2)),
    mcc = map_dbl(mat, ~mcc(confusionM = .x)),
    acc = map_dbl(mat, accuracy)
  )

testing_mcc2 <- tibble(n_righty_correct = seq(0, 200, 10)) %>%
  mutate(
    mat = map(n_righty_correct,
              ~matrix(data = c(20, 200 - .x, 0, .x), nrow = 2)),
    mcc = map_dbl(mat, ~mcc(confusionM = .x)),
    acc = map_dbl(mat, accuracy)
  )

testing_mcc3 <- expand_grid(n_lefty_correct = seq(0, 20, 2),
                            n_righty_correct = seq(0, 200, 10)) %>%
  mutate(
    mat = map2(n_lefty_correct, n_righty_correct,
              ~matrix(data = c(.x, 20 - .x, 200 - .y, .y), nrow = 2)),
    mcc = map_dbl(mat, ~mcc(confusionM = .x)),
    acc = map_dbl(mat, accuracy),

    lefty8090 = case_match(n_lefty_correct,
                          c(16, 18) ~ n_lefty_correct,
                          TRUE ~ 0)

  )

testing_values <- bind_rows(testing_mcc, testing_mcc2) %>%
  select(-mat) %>%
  pivot_longer(c(mcc, acc)) %>%
  mutate(
    x = if_else(is.na(n_righty_correct), n_lefty_correct / 20,
                n_righty_correct / 200),
    group = if_else(is.na(n_righty_correct), "L", "R")
  )

png("plots/mcc.png", width = 6, height = 4, units = "in", res = 300)

ggplot(testing_values, aes(x = x, y = value, color =  group)) +
  geom_point(aes(shape = name)) +
  geom_line(aes(linetype = name)) +
  theme_bw() +
  labs(x = "Correct (%)", y = "Value", color = "Metric")

dev.off()

ggplot(testing_mcc3, aes(x = n_righty_correct / 200 * 100, y = mcc,
                         color = as.factor(lefty8090))) +
  geom_point(size = 1) +
  geom_line(aes(group = n_lefty_correct), linewidth = 0.25) +
  geom_vline(xintercept = c(80, 90), color = c("orange", "darkred")) +
  scale_x_continuous(breaks = seq(0, 100, 20)) +
  scale_y_continuous(breaks = seq(-1, 1, by = 0.2)) +
  scale_color_manual(values = c("orange", "darkred"), na.value = "grey",
                     labels = c("80%", "90%")) +
  theme_bw() +
  labs(x = "Dextral correct %", color = "Sinistral correct %",
       y = "MCC")

testing_mcc3 %>%
  filter(
    n_righty_correct %in% (200 * c(0.5, 0.8, 0.9)),
    n_lefty_correct %in% (20 * c(0.5, 0.8, 0.9)),
  )
