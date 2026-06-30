library(tidyverse)

hemi <- read_csv("/Users/tkmd/MyDrive/Projects/hemisphere_fingerprinting/code/modeling_new/sizes_mcc.csv",
                     show_col_types = FALSE)

gender <- read_csv("/Users/tkmd/MyDrive/Projects/hemisphere_fingerprinting/code/modeling_new/sizes_mcc_gender.csv",
                     show_col_types = FALSE) 

age <- read_csv("/Users/tkmd/MyDrive/Projects/hemisphere_fingerprinting/code/modeling_new/sizes_mcc_age.csv",
                     show_col_types = FALSE)

hand <- read_csv("/Users/tkmd/MyDrive/Projects/hemisphere_fingerprinting/code/modeling_new/sizes_mcc_hand.csv",
                     show_col_types = FALSE)

all_results <- bind_rows(hemi, gender, age, hand)


ggplot(all_results, aes(x = size)) +
  geom_ribbon(aes(ymin = mcc_ci_lo, ymax = mcc_ci_hi, fill = outcome), alpha = 0.5) +
  geom_line(aes(y = mcc, color = outcome), linewidth = 1) +
  scale_x_continuous(breaks = results$size, transform = "log2") +
  theme_bw() +
  labs(title = "MCC by outcome, log x", x = "Size", y = "MCC")

ggplot(all_results, aes(x = size)) +
  geom_ribbon(aes(ymin = mcc_ci_lo, ymax = mcc_ci_hi, fill = outcome), alpha = 0.5) +
  geom_line(aes(y = mcc, color = outcome), linewidth = 1) +
  theme_bw() +
  labs(title = "MCC by outcome, linear x", x = "Size", y = "MCC")
