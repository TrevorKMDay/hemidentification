library(tidyverse)
library(mltools)
library(doParallel)
library(foreach)
library(rstatix)

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

hc <- qs::qread("../inputs/hemiconnectome.rds")

demo <- hc %>%
  select(sub, handedness) %>%
  distinct()

table(demo$handedness)

all_files1 <- list.files("../results/", pattern = ".*_data-.*.csv",
                         full.names = TRUE, recursive = TRUE)

lefty_files <- str_subset(all_files1, "_hands-lefty_") %>%
  str_subset("osampl", negate = TRUE)

hand_files <- str_subset(all_files1, "outcome-hand")
subset_files <- str_subset(all_files1, "hands-righty_data-base_subset-74")

results <- tibble(f = c(lefty_files, hand_files, subset_files)) %>%
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
  select(-`...1`, -f)

results_nested <- results %>%
  mutate(

    n = map_int(data, nrow),

    outcome = str_extract(f, "outcome-[^_]*"),
    method = str_extract(f, "method-[^_]*"),
    group = str_extract(f, "test-[^_]*"),
    input = str_extract(f, "data-[^_]*"),

    subset = str_extract(f, "subset-[^_]*") %>%
      replace_na("none"),

    across(c(outcome, method, group, input, subset), ~str_remove(.x, "^.*-"))

  ) %>%
  unnest(data) %>%
  mutate(
    ground_hand = str_remove(ground, ".H$"),
    ground_hemi = str_extract(ground, ".H$"),
    pred_hand = str_remove(predicted, ".H$"),
    pred_hemi = str_extract(predicted, ".H$"),
  ) %>%
  group_by(outcome, method, group, input, subset, n) %>%
  select(-`...1`, -f) %>%
  mutate(
    input = if_else(subset != "none", paste0(input, subset), input)
  )

subset <- results_nested %>%
  filter(
    input == "base74"
  )

results_mcc <- results_nested %>%
  nest() %>%
  mutate(

    hand = map(data, bootstrap_mcc,
                R = 1000, ground = "ground_hand", predicted = "pred_hand",
                levels = c("righty", "lefty"), method = 2,
              .progress = TRUE),

    # Specify only evaluating hemi
    hemi = map(data, bootstrap_mcc,
                   R = 1000, ground = "ground_hemi", predicted = "pred_hemi",
                   levels = c("LH", "RH"), method = 2,
                   .progress = TRUE),

    acc = map_dbl(data, ~sum(.$ground_hemi == .$pred_hemi) / nrow(.)),

    n_errors = map_int(data, ~sum(.$ground_hemi != .$pred_hemi))

  ) %>%
  unnest(c(hand, hemi), names_sep = ".")

# results_mcc_summary <- results_mcc %>%
#   ungroup() %>%
#   group_by(method, input) %>%
#   summarize(
#     mcc_min = min(mcc_mean),
#     mcc_max = max(mcc_mean),
#     acc_min = min(acc),
#     acc_max = max(acc),
#   ) %>%
#   mutate(
#     across(where(is.numeric), ~round(., 3))
#   )

# Compare lefty and righty results ====

righty_results <- read_rds("right_handed_results.rds") %>%
  ungroup() %>%
  select(method, group, input, subset, mcc_mean) %>%
  rename(
    righty_mcc = mcc_mean
  )

lefty_v_righty <- results_mcc %>%
  ungroup() %>%
  filter(
    input %in% c("base", "langtask"),
    outcome == "hemi"
  ) %>%
  select(method, group, input, subset, hemi.mcc_mean) %>%
  rename(
    lefty_mcc = hemi.mcc_mean
  ) %>%
  left_join(righty_results, by = join_by(method, group, input, subset))

lefty_v_righty_test <- lefty_v_righty %>%
  group_by(method, input) %>%
  nest() %>%
  mutate(
    n = map_int(data, nrow),
    mean_lefty = map_dbl(data, ~mean(.x$lefty_mcc)),
    mean_righty = map_dbl(data, ~mean(.x$righty_mcc))
  )

lefty_v_righty_test2 <- lefty_v_righty_test %>%
  filter(
    !(mean_lefty == 1 &  mean_righty == 1)
  ) %>%
  mutate(

    t_test = if_else(mean_righty == 1,
                     map(data, ~t.test(.x$lefty_mcc, mu = 1)),
                     map(data, ~t.test(.x$lefty_mcc, .x$righty_mcc))),

    p_value = map_dbl(t_test, ~ .x$p.value)

  )

lefty_v_righty_test2$p_cor = p.adjust(lefty_v_righty_test2$p_value,
                                      method = "fdr", n = 6)

# Compare small group ====

lefties_base <- results_mcc %>%
  ungroup() %>%
  filter(
    outcome == "hemi",
    input == "base"
  ) %>%
  select(method, group, hemi.mcc_mean) %>%
  rename(
    all_lefties = hemi.mcc_mean
  )

righties_base <- righty_results %>%
  filter(
    input == "base"
  )

hc92 <- results_mcc %>%
  ungroup() %>%
  filter(
    input == "base92"
  ) %>%
  select(method, group, hemi.mcc_mean) %>%
  rename(
    hc92 = hemi.mcc_mean
  ) %>%
  left_join(lefties_base) %>%
  left_join(righties_base)

t.test(hc92$hc92, hc92$righty_mcc, paired = "TRUE")

hc92_models <- hc92 %>%
  group_by(method) %>%
  nest() %>%
  mutate(
    n = map_int(data, nrow),
    mean_lefties = map_dbl(data, ~mean(.x$all_lefties)),
    mean_hc92 = map_dbl(data, ~mean(.x$hc92)),

    t_test = if_else(mean_lefties == 1,
                     map(data, ~t.test(.x$hc92, mu = 1)),
                     map(data, ~t.test(.x$all_lefties, .x$hc92))),

    p_value = map_dbl(t_test, ~ .x$p.value)

  )

hc92_models$p_adj <- p.adjust(hc92_models$p_value, "fdr", n = 3)

# Hemi lefties =====

hemi_results <- results_mcc %>%
  ungroup() %>%
  filter(
    outcome == "hemi" | input == "base74",
    input != "full"
  ) %>%
  select(-starts_with("hands"))

ggplot(hemi_results, aes(x = input, y = hemi.mcc_mean, color = group)) +
  geom_hline(yintercept = c(0.4, 0.6), linetype = "dashed", color = "red") +
  geom_pointrange(aes(ymin = hemi.mcc_ci_lo, ymax = hemi.mcc_ci_hi),
                  position = position_dodge2(width = 0.5),
                  alpha = 0.75) +
  scale_x_discrete(labels = c("HC\n(L)", "HC74\n(R)", "LgHC\n(L)",
                              "HC-R\n(L)")) +
  scale_y_continuous(limits = c(NA, NA), breaks = c(0.5, 0.75, 1)) +
  facet_wrap(vars(method)) +
  theme_bw() +
  labs(x = "Input data", y = "MCC (95% CI)")

ggsave("plots/lefties_MCC.png", width = 6.5, height = 3.5, units = "in")

# hemi_aov1 <- aov(hemi.mcc_mean ~ input + method + Error(group/input)
#                  data = hemi_results)
#
# summary(hemi_aov1)

# TukeyHSD(hemi_aov1, which = c("input", "method"))

hemi_aov1 <- anova_test(hemi.mcc_mean ~ input + method + group, wid = group,
                         between = c(input, method),
                         data = hemi_results)

hemi_aov1

# Handedness classification

results_mcc %>%
  filter(
    input != "base74"
  ) %>%
  pull(hand.mcc_mean) %>%
  range()

results_mcc %>%
  filter(
    input %in% c("base", "langtask", "osampl")
  ) %>%
  pull(hemi.mcc_mean) %>%
  range()

hand_results <- results_mcc %>%
  filter(
    outcome == "hand",
    !(input %in% c("base74", "rev"))
  )


ggplot(hand_results, aes(x = input, y = hand.mcc_mean, color = group)) +
  geom_hline(yintercept = 0, linetype = "solid", color = "red") +
  geom_hline(yintercept = 0.4, linetype = "dashed", color = "red") +
  geom_pointrange(aes(ymin = hand.mcc_ci_lo, ymax = hand.mcc_ci_hi),
                  position = position_dodge2(width = 0.5),
                  alpha = 0.75) +
  scale_x_discrete(labels = c("HC\n(L)", "FC\n(L)", "LgHC\n(L)", "HCO\n(L)",
                              "TC\n(L)")) +
  facet_wrap(vars(method)) +
  theme_bw() +
  labs(x = "Input data", y = "MCC (95% CI)")

results_mcc %>%
  filter(
    input == "base92"
  ) %>%
  group_by(method) %>%
  summarize(
    min = min(hemi.mcc_mean),
    max = max(hemi.mcc_mean)
  )

ggsave("plots/handedness_MCC.png", r_hands, width = 6.5, height = 3.5,
       units = "in")

# ANOVA ====

hand_aov1 <- anova_test(hand.mcc_mean ~ input + method + group, wid = group,
                        between = c(input, method),
                        data = ungroup(hand_results))

hand_aov1

# results_anova <- aov(hand.mcc_mean ~ method + group + input,
#                      data = results_nosubsample)
#
# TukeyHSD(results_anova, which = c("group", "input"))

# LDA results ====

LDA <- read_csv("../lda_full/lda_class_xform.csv", show_col_types = FALSE) %>%
  select(sub, hemi, handedness, EHI, `0`, `1`, `2`) %>%
  rename(
    LD1 = `0`,
    LD2 = `1`,
    LD3 = `2`
  )

LDA_wide <- LDA %>%
  pivot_wider(names_from = hemi, values_from = starts_with("LD"))

cors1 <- LDA_wide %>%
  group_by(handedness) %>%
  nest() %>%
  mutate(
    rcorr = map(data, ~Hmisc::rcorr(as.matrix(.x[, c("LD1_RH", "LD1_LH")]))),
    r = map_dbl(rcorr, ~.x$r[1, 2]),
    p = map_dbl(rcorr, ~.x$P[1, 2])
  )

