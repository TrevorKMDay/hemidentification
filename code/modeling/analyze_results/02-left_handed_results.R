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

setwd(paste0("~/MyDrive/Projects/hemisphere_fingerprinting/code/",
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

results_file <- "results_mcc-lefties.rds"

if (!file.exists(results_file)) {

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

  qs::qsave(results_mcc, results_file)

} else {

  results_mcc <- qs::qread(results_file)

}


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
    input %in% c("base", "langtask")
  ) %>%
  select(method, group, input, hemi.mcc_mean) %>%
  pivot_wider(values_from = hemi.mcc_mean, names_from = input,
              names_prefix = "lefties.")

righties_base <- righty_results %>%
  filter(
    input %in% c("base", "langtask")
  ) %>%
  pivot_wider(values_from = righty_mcc, names_from = input,
              names_prefix = "righty_mcc.")

hc_small <- results_mcc %>%
  ungroup() %>%
  filter(
    input == "base74"
  ) %>%
  select(method, group, hemi.mcc_mean) %>%
  rename(
    hc_small.base = hemi.mcc_mean
  ) %>%
  left_join(lefties_base, by = join_by(method, group)) %>%
  left_join(righties_base, by = join_by(method, group))

## Small/big R =====

t.test(hc_small$hc_small.base, hc_small$righty_mcc.base, paired = TRUE)

subset_d <- (mean(hc_small$righty_mcc) - mean(hc_small$hc_small)) /
              sd(c(hc_small$righty_mcc, hc_small$hc_small))

# data:  hc_small$hc_small and hc_small$righty_mcc
# t = -4.2188, df = 14, p-value = 0.0008587
# alternative hypothesis: true mean difference is not equal to 0
# 95 percent confidence interval:
#   -0.018131909 -0.005909451
# sample estimates:
#   mean difference
# -0.01202068

hc_small_models <- hc_small %>%
  select(method, group, ends_with("base"), ends_with("langtask"),
         -hc_small.base) %>%
  pivot_longer(-c(method, group), names_to = c("hand", "input"),
               names_sep = "[.]", values_to = "mcc") %>%
  mutate(
    group2 = if_else(hand == "lefties", paste0(group, "L"), paste0(group, "R"))
  )

lVr_lm <- lmerTest::lmer(mcc ~ hand*method*input + (1|group2),
                           data = hc_small_models)

summary(lVr_lm)

## HC-R ====

lefties_rev <- results_mcc %>%
  ungroup() %>%
  filter(
    input == "rev"
  ) %>%
  select(method, group, input, hemi.mcc_mean)

rev_lm <- lmerTest::lmer(hemi.mcc_mean ~ method + (1|group),
                         data = lefties_rev)

summary(rev_lm)

# Hemi lefties =====

hemi_results <- results_mcc %>%
  ungroup() %>%
  filter(
    outcome == "hemi" | input == "base74",
    input != "full",
    input != "baseZ"
  ) %>%
  select(-starts_with("hands"))

ggplot(hemi_results, aes(x = input, y = hemi.mcc_mean, color = group)) +
  geom_hline(yintercept = c(0.4, 0.6), linetype = "dashed", color = "red") +
  geom_pointrange(aes(ymin = hemi.mcc_ci_lo, ymax = hemi.mcc_ci_hi),
                  position = position_dodge2(width = 0.5),
                  alpha = 0.75) +
  facet_wrap(vars(method)) +
  scale_x_discrete(labels = c("HC\n(L)", "HC74\n(R)", "LgHC\n(L)",
                              "HC-R\n(L)")) +
  scale_y_continuous(limits = c(NA, NA), breaks = c(0.5, 0.75, 1)) +
  theme_bw() +
  labs(x = "Input data", y = "MCC (95% CI)", color = "Test fold")

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
  labs(x = "Input data", y = "MCC (95% CI)", color = "Test fold")

results_mcc %>%
  filter(
    input == "base92"
  ) %>%
  group_by(method) %>%
  summarize(
    min = min(hemi.mcc_mean),
    max = max(hemi.mcc_mean)
  )

ggsave("plots/handedness_MCC.png", width = 6.5, height = 3.5,
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

