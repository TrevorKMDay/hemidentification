library(tidyverse)
library(mltools)
library(doParallel)
library(foreach)
library(patchwork)

select <- dplyr::select

# Parallel setup ====

cores <- detectCores()
cl <- makeCluster(cores[1] - 1) # not to overload your computer
registerDoParallel(cl)

# Setup ====

setwd(paste0("~/MyDrive/Projects/hemisphere_fingerprinting/code/modeling/",
              "analyze_results/"))

source("bootstrap_mcc.R")

# There are a lot of base variations, including different hand training splits
# (righty2/nonrighty/hand2), subsets for comparision (subset), and predicting
# handedness, not hemisphere (outcome-hand); remove those

all_files1 <- list.files("../results/", pattern = ".*_data-.*.csv",
                         full.names = TRUE, recursive = TRUE) %>%
  str_subset("/(base|connNormal)/") %>%
  str_subset("nonrighty|righty2|hand2|subset|outcome-hand", negate = TRUE)

base_files <- str_subset(all_files1, "base")

# Read data in from files
results <- tibble(f = all_files1) %>%
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

    across(c(method, group, input), ~str_remove(.x, "^.*-"))

  ) %>%
  unnest(data) %>%
  mutate(
    ground_hemi = str_extract(ground, ".H$"),
    pred_hemi = str_extract(predicted, ".H$"),
  ) %>%
  group_by(method, group, input, n) %>%
  filter(
    ground %in% c("rightyRH", "rightyLH", "RH", "LH")
  ) %>%
  select(-`...1`, -f)

results_nested <- results %>%
  mutate(

    n = map_int(data, nrow),

    method = str_extract(f, "method-[^_]*"),
    group = str_extract(f, "test-[^_]*"),
    input = str_extract(f, "data-[^_]*"),
    hands = str_extract(f, "hands-[^_]*"),
    input = str_extract(f, "data-[^_]*"),

    across(c(method, group, input, hands, input), ~str_remove(.x, "^.*-"))

  ) %>%
  unnest(data) %>%
  mutate(
    ground_hemi = str_extract(ground, ".H$"),
    pred_hemi = str_extract(predicted, ".H$"),
  ) %>%
  group_by(method, group, hands, input, n) %>%
  filter(
    # ground %in% c("rightyRH", "rightyLH", "rightyRH", "rightyLH", "LH", "RH"),
  ) %>%
  select(-`...1`, -f)

results_mcc <- results_nested %>%
    nest() %>%
    mutate(

        # Specify only evaluating hemi
        mcc = map(data, bootstrap_mcc,
                    R = 1000, ground = "ground_hemi", predicted = "pred_hemi",
                    levels = c("LH", "RH"), method = 2,
                .progress = TRUE),

        acc = map_dbl(data, ~sum(.$ground_hemi == .$pred_hemi) / nrow(.)),

        n_errors = map_int(data, ~sum(.$ground_hemi != .$pred_hemi))

    ) %>%
    unnest(mcc)

ggplot(results_mcc, aes(x = interaction(input, method), y = mcc_mean, color = group)) +
  geom_pointrange(aes(ymin = mcc_ci_lo, ymax = mcc_ci_hi, shape = input),
                  position = position_dodge2(width = 0.5),
                  alpha = 0.75) +
  scale_y_continuous(limits = c(0.5, NA)) +
  facet_wrap(vars(hands)) +
  theme_bw() +
  labs(x = "Method", y = "MCC", color = "Test\ngroup") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

results_wide <- results_mcc %>%
  ungroup() %>%
  select(method, group, input, hands, mcc_mean) %>%
  pivot_wider(names_from = input, values_from = mcc_mean)

t.test(results_wide$base, results_wide$connNormal, paired = TRUE)


# t = 1.8994, df = 44, p-value = 0.06408
# alternative hypothesis: true mean difference is not equal to 0
# 95 percent confidence interval:
#  -0.0004131544  0.0139499308
# sample estimates:
# mean difference
#     0.006768388

# base > connNormal, at p = .064