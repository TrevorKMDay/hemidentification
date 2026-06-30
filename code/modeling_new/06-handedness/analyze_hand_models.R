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
  str_subset("gender", negate = TRUE)

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
  unnest(data) 

outcome_hand <- data %>%
  left_join(folds) %>%
  group_by(clf, input, outcome, group) %>%
  nest() %>%
  filter(
    outcome == "outcome-handedness"
  ) %>%
  mutate(
    # Bootstrap MCC
    mcc_outcome = "mcc_hand",
    mcc_hand = map(data, bootstrap_mcc, ground = "handedness")
  ) %>%
  unnest(mcc_hand)

outcome_hand_summary <- data %>%
  left_join(folds) %>%
  group_by(clf, input, outcome) %>%
  nest() %>%
  filter(
    outcome == "outcome-handedness"
  ) %>%
  mutate(
    # Bootstrap MCC
    mcc_outcome = "mcc_hand",
    mcc_hand = map(data, bootstrap_mcc, ground = "handedness")
  ) %>%
  unnest(mcc_hand)

class1 <- data %>%
  left_join(folds) %>%
  filter(
    outcome == "outcome-class"
  ) %>%
  mutate(

    hemi_pred = str_extract(predicted, "[LR]H$"),
    hand_pred = str_remove(predicted, "[LR]H$")

  ) 

outcome_class <- class1 %>%
  group_by(clf, input, outcome, group)  %>%
  nest() %>%
  mutate(

    # Bootstrap MCC
    mcc_overall = map(data, bootstrap_mcc, ground = "class",
                      .progress = TRUE),

    # Bootstrap MCC
    mcc_hemi = map(data, bootstrap_mcc, predicted = "hemi_pred", 
                        ground = "hemi",
                      .progress = TRUE),

    mcc_hand = map(data, bootstrap_mcc, predicted = "hand_pred", 
                        ground = "handedness",
                      .progress = TRUE),


  ) %>%
  unnest(starts_with("mcc_"), names_sep = ".")

outcome_class2 <- outcome_class %>% 
  pivot_longer(starts_with("mcc"), names_sep = "[.]", 
               names_to = c("mcc_outcome", "msmt")) %>%
  pivot_wider(names_from = msmt)

outcome_class_summary <- class1 %>%
  group_by(clf, input, outcome)  %>%
  nest() %>%
  mutate(

    # Bootstrap MCC
    mcc_overall = map(data, bootstrap_mcc, ground = "class",
                      .progress = TRUE),

    # Bootstrap MCC
    mcc_hemi = map(data, bootstrap_mcc, predicted = "hemi_pred", 
                        ground = "hemi",
                      .progress = TRUE),

    mcc_hand = map(data, bootstrap_mcc, predicted = "hand_pred", 
                        ground = "handedness",
                      .progress = TRUE),


  ) %>%
  unnest(starts_with("mcc_"), names_sep = ".")

outcome_class_summary2 <- outcome_class_summary %>% 
  pivot_longer(starts_with("mcc"), names_sep = "[.]", 
               names_to = c("mcc_outcome", "msmt")) %>%
  pivot_wider(names_from = msmt)

# Put together for plotting ====

all_outcome <- bind_rows(outcome_hand, outcome_class2) %>%
  mutate(

    input = str_remove(input, "data-") %>%
      recode_values(
        "fc" ~ "FC",
        "trans" ~ "TC",
        "hc" ~ "HC",
        "hco" ~ "HCO",
        "lang" ~ "LgHC"
      ),

    mcc_outcome = mcc_outcome %>%
      str_remove("mcc_") %>%
      recode_values(
        "hand" ~ "Handedness",
        "hemi" ~ "Hemisphere",
        "overall" ~ "Overall"
      )

  )

all_summaries <- bind_rows(outcome_hand_summary, outcome_class_summary2) %>%
  mutate(

    input = str_remove(input, "data-") %>%
      recode_values(
        "fc" ~ "FC",
        "trans" ~ "TC",
        "hc" ~ "HC",
        "hco" ~ "HCO",
        "lang" ~ "LgHC"
      ),

    mcc_outcome = mcc_outcome %>%
      str_remove("mcc_") %>%
      recode_values(
        "hand" ~ "Handedness",
        "hemi" ~ "Hemisphere",
        "overall" ~ "Overall"
      )

  )

all_lda_outcomes <- filter(all_outcome, clf == "clf-lda")
all_lda_summaries <- filter(all_summaries, clf == "clf-lda")

ggplot(all_lda_outcomes, aes(x = input, color = mcc_outcome)) +
  geom_pointrange(
    aes(y = mcc_mean, ymin = mcc_ci_lo, ymax = mcc_ci_hi),
    position = position_dodge2(width = 0.25),
    alpha = 0.75
  ) +
  geom_pointrange(
    data = all_lda_summaries,
    aes(x = input, y = mcc_mean, ymin = mcc_ci_lo, ymax = mcc_ci_hi, 
        fill = mcc_outcome),
    position = position_dodge2(width = 0.25),
    size = 0.9, shape = 21, color = "black"
  ) +
  geom_hline(yintercept = 0, color = "red") +
  scale_x_discrete(limits = c("HC", "HCO", "LgHC", "FC", "TC")) +
  theme_bw() +
  labs(x = "Input dataset", y = "MCC (95% CI)", color = "Outcome", 
        fill = "Outcome") +
  theme(legend.position = "bottom")

ggsave("plots/fourway_models_mcc.png", width = 6.5, height = 4.5)

other_outcomes <- all_outcome %>%
  filter(
    clf != "clf-lda"
  ) %>%
  mutate(
    clf = toupper(str_remove(clf, "clf-"))
  )


other_summaries <- all_summaries %>%
  filter(
    clf != "clf-lda"
  ) %>%
  mutate(
    clf = toupper(str_remove(clf, "clf-"))
  )


ggplot(other_outcomes, aes(x = input, color = mcc_outcome)) +
  geom_pointrange(
    aes(y = mcc_mean, ymin = mcc_ci_lo, ymax = mcc_ci_hi),
    position = position_dodge2(width = 0.25),
    alpha = 0.75
  ) +
  geom_pointrange(
    data = other_summaries,
    aes(x = input, y = mcc_mean, ymin = mcc_ci_lo, ymax = mcc_ci_hi, 
        fill = mcc_outcome),
    position = position_dodge2(width = 0.25),
    size = 0.9, shape = 21, color = "black"
  ) +
  geom_hline(yintercept = 0, color = "red") +
  scale_x_discrete(limits = c("HC", "HCO", "LgHC", "FC", "TC")) +
  facet_wrap(vars(clf)) +
  theme_bw() +
  labs(x = "Input dataset", y = "MCC (95% CI)", color = "Outcome", 
        fill = "Outcome") +
  theme(legend.position = "bottom")

ggsave("plots/fourway_models_NNandSVC_mcc.png", width = 6.5, height = 4.5)

# summary values

range(all_lda_summaries$mcc_mean[all_lda_summaries$mcc_outcome == "Handedness"])

range(all_summaries$mcc_mean[all_summaries$mcc_outcome == "Handedness" &
                              all_summaries$clf == "clf-nn"])

range(all_summaries$mcc_mean[all_summaries$mcc_outcome == "Handedness" &
                              all_summaries$clf == "clf-svc"])


range(all_summaries$mcc_mean[all_summaries$mcc_outcome == "Hemisphere" &
                              all_summaries$clf == "clf-nn"])

range(all_summaries$mcc_mean[all_summaries$mcc_outcome == "Hemisphere" &
                              all_summaries$clf == "clf-svc"])

all_hands <- all_lda_outcomes %>%
  filter(
    mcc_outcome == "Handedness"
  )

hand_anova <- aov(mcc_mean ~ group + input, data = all_hands)

hand_tukey <- TukeyHSD(hand_anova)

hand_tukey$group %>%
  as_tibble(rownames = "comp") %>%
  filter(
    `p adj` < .05
  )

table(all_hands$mcc_mean > 0.4)
