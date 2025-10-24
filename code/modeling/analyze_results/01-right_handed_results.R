library(tidyverse)
library(mltools)
library(doParallel)
library(foreach)

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

all_files1 <- list.files("../results/", pattern = ".*_data-.*.csv",
                         full.names = TRUE, recursive = TRUE)

righty_files <- str_subset(all_files1, "_hands-righty_") %>%
  str_subset("osampl", negate = TRUE)

# Read data in from files
results <- tibble(f = righty_files) %>%
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

    subset = str_extract(f, "subset-[^_]*") %>%
      replace_na("none"),

    across(c(method, group, input, subset), ~str_remove(.x, "^.*-"))

  ) %>%
  unnest(data) %>%
  mutate(
    ground_hemi = str_extract(ground, ".H$"),
    pred_hemi = str_extract(predicted, ".H$"),
  ) %>%
  group_by(method, group, input, subset, n) %>%
  filter(
    ground %in% c("rightyRH", "rightyLH", "LH", "RH"),
    subset == "none"
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

results_mcc_summary <- results_mcc %>%
  ungroup() %>%
  group_by(method, input) %>%
  summarize(
    mcc_min = min(mcc_mean),
    mcc_max = max(mcc_mean),
    acc_min = min(acc),
    acc_max = max(acc),
  ) %>%
  mutate(
    across(where(is.numeric), ~round(., 3))
  )

ggplot(results_mcc, aes(x = input, y = mcc_mean, color = group)) +
  geom_pointrange(aes(ymin = mcc_ci_lo, ymax = mcc_ci_hi),
                  position = position_dodge2(width = 0.5),
                  alpha = 0.75) +
  scale_y_continuous(limits = c(0.9, 1.01),
                     breaks = c(0.9, 0.95, 1.0)) +
  scale_x_discrete(labels = c("HC", "LgHC", "HC-R")) +
  facet_wrap(vars(method)) +
  theme_bw() +
  labs(x = "Input data", y = "MCC (95% CI)")

ggsave("plots/righties_MCC.png", width = 5, height = 3.5, units = "in")

write_rds(results_mcc, "right_handed_results.rds")

# LDA ====

righty_LDA <- read_csv("../righty_full.csv", show_col_types = FALSE) %>%
  select(sub, gender, hemi, predicted, `0`) %>%
  rename(
    LD1 = `0`
  )

righty_LDA_wide <- righty_LDA %>%
  pivot_wider(id_cols = c(sub, gender), names_from = hemi,
              values_from = LD1) %>%
  mutate(
    diff = RH - LH
  )

ggplot(righty_LDA_wide, aes(x = diff, fill = gender)) +
  geom_density()

Hmisc::rcorr(as.matrix(righty_LDA_wide[, c("RH", "LH")]))

t.test(righty_LDA_wide$diff[righty_LDA_wide$gender == "F"],
       righty_LDA_wide$diff[righty_LDA_wide$gender == "M"])

ggplot(righty_LDA_wide, aes(x = LH, y = RH)) +
  geom_point()

righty_LDA_plot <- ggplot(righty_LDA, aes(x = LD1, y = sub, color = hemi)) +
  geom_point(alpha = 0.5) +
  scale_y_discrete(labels = NULL, breaks = NULL,
                   expand = expansion(mult = 0.05)) +
  theme_bw() +
  labs(y = "Participant", color = "Hemisphere") +
  theme(legend.position = "bottom")

snl_plot <- righty_LDA_plot +
  theme(text = element_text(size = 24))

ggsave("plots/snl_righty_LD_plot.png", snl_plot, width = 10, height = 8,
       units = "in")

# Scalings ====

hc <- qs::qread("../inputs/hemiconnectome.rds")
enrich <- read_rds("../../enrichment/hemisphere_model_ps.rds")

hc_long <- hc %>%
  select(sub, hemi, contains("_"), -age_group) %>%
  pivot_longer(-c(sub, hemi), names_to = "feature") %>%
  pivot_wider(names_from = hemi) %>%
  mutate(
    diff = LH -RH
  )

hc_avg_diff <- hc_long %>%
  group_by(feature) %>%
  summarize(
    n = n(),
    avg_diff = mean(diff)
  )

scalings <- read_csv("../lda_full/full_scalings.csv",
                     show_col_types = FALSE) %>%
  left_join(hc_avg_diff, by = join_by(feature)) %>%
  left_join(enrich, by = join_by(feature, LD1)) %>%
  mutate(
    sig = p < .01,
    greater = if_else(avg_diff > 0, "LH", "RH")
  )

table(scalings$sig, scalings$greater)


Hmisc::rcorr(as.matrix(scalings[, c("LD1", "avg_diff")]))

ggplot(scalings, aes(x = LD1, y = avg_diff, color = p < .01)) +
  geom_point(alpha = 0.1)
