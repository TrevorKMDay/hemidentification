library(tidyverse)
library(mltools)

# Example confusion matrix  ====

# Use "confusion matrix" instead of TP/FP, etc. for clarity

# actual      LH RH
# assigned LH  .  .
#          RH  .  .

h1 <- matrix(data = c(45, 5, 4, 46), nrow = 2)
mcc(confusionM = h1)

gen_h1_confmatrix <- function(LH_pct, RH_pct, n = 100, jitter = 0.2) {

  LH_pct2 <- (LH_pct + runif(1, -jitter, jitter))
  LH_pct2 <- if_else(LH_pct2 > 1, 1, if_else(LH_pct2 < 0, 0, LH_pct2))

  RH_pct2 <- (RH_pct + runif(1, -jitter, jitter))
  RH_pct2 <- if_else(RH_pct2 > 1, 1, if_else(RH_pct2 < 0, 0, RH_pct2))

  LH1 <- round(100 * LH_pct2)
  LH2 <- 100 - LH1

  RH1 <- round(100 * RH_pct2)
  RH2 <- 100 - RH1

  mat <- matrix(data = c(LH1, LH2, RH2, RH1), nrow = 2)

  return(mat)

}

repeats <- expand_grid(fold = LETTERS[1:5],
                       tasks = c("ALL8", "RS1"),
                       model = c("LDA", "SVM", "NN"),
                       hemis = c("linked", "one")) %>%
  rowwise() %>%
  mutate(
    LH = case_match(fold,
                    "A" ~ 0.95, "B" ~ 0.9, "C" ~ 0.85, "D" ~ 0.8, "E" ~ 0.75),

    RH = LH - runif(1, -0.025, 0.025),
  ) %>%
  ungroup() %>%
  mutate(
    data = map2(LH, RH, ~gen_h1_confmatrix(.x, .y)),
    mcc = map_dbl(data, ~mcc(confusionM = .x))
  )

ggplot(repeats, aes(x = model, y = mcc, color = model, shape = tasks)) +
  geom_point(position = position_dodge(width = 0.25)) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw()

lm(mcc ~ tasks + model + hemis + fold, data = repeats) %>%
  summary()

library(lme4)

lmer(mcc ~ tasks + model + hemis + (1|fold), data = repeats)

mcc(confusionM = matrix(c(45, 4, 4, 4, 45, 4, 4, 4, 45), ncol = 3))

# Macroaverage MCC ====

big_matrix <- matrix(
    c(10, 1, 1, 1,
      1, 11, 1, 1,
      1, 1, 12, 1,
      1, 1, 1, 13),
    ncol = 4, byrow = TRUE
  )

ma_mcc <- function(conf_matrix, ij) {

  conf_diag <- diag(conf_matrix)
  TP <- conf_diag[ij]
  TO <- sum(conf_diag[-ij])

  the_row <- conf_matrix[ij, ]
  FP <- sum(the_row[-ij])

  the_col <- conf_matrix[, ij]
  FN <- sum(the_col[-ij])

  mcc_val <- mcc(TP = TP, FP = FP, TN = TO, FN = FN)

  return(mcc_val)

}

ma_mcc(big_matrix, 1)
mcc(confusionM = big_matrix)

mcc(confusionM = big_matrix[1:2, 1:2])

rainio_table2 <- matrix(
    c(120, 7, 9, 4,
      15, 116, 3, 6,
      12, 13, 115, 0,
      2, 96, 4, 38),
    ncol = 4, byrow = TRUE
  )

mcc(confusionM = rainio_table2)
