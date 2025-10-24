library(tidyverse)
library(mltools)
library(doParallel)
library(foreach)

select <- dplyr::select

# Parallel setup ====

cores <- detectCores()
cl <- makeCluster(cores[1] - 1) #not to overload your computer
registerDoParallel(cl)


# Setup ====

setwd("~/Google Drive/My Drive/Projects/hemisphere_fingerprinting/code/modeling/analyze_results/")
source("bootstrap_mcc.R")

# Find data ====

all_files1 <- list.files("../results/", pattern = ".*_data-.*.csv",
                         full.names = TRUE, recursive = TRUE)

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
    outcome = str_extract(f, "outcome-[^_]*"),
    group = str_extract(f, "test-[^_]*"),
    input = str_extract(f, "data-[^_]*"),

    subset = str_extract(f, "subset-[^_]*"),

    hands = replace_na(str_extract(f, "hands-[^_]*"), "all"),
    hemis = replace_na(str_extract(f, "hemi-[^_]*"), "all"),

    across(c(outcome, method, group, hands, hemis, input, subset),
           ~str_remove(.x, "^.*-"))

  ) %>%
  unnest(data) %>%
  mutate(
    ground_hand = str_remove(ground, ".H$"),
    ground_hemi = str_extract(ground, ".H$"),
    pred_hand = str_remove(predicted, ".H$"),
    pred_hemi = str_extract(predicted, ".H$"),
  ) %>%
  select(-`...1`, -f) %>%
  group_by(outcome, method, group, hands, hemis, input, subset, n)

results_big %>%
  filter(
    !(ground_hemi %in% c("LH", "RH"))
  )

results2 <- results_big %>%
  nest() %>%
  mutate(
    input = if_else(!is.na(subset), paste0(input, subset), input)
  ) %>%
  filter(
    input != "base24"
  )

# Results by outcome ====

## Hemisphere ====

results_hemis <- results2 %>%
  filter(
    hemis == "all",
    !(input %in% c("trans", "full")),
    outcome == "hemi"
  ) %>%
  mutate(

    # Specify only evaluating hemi
    mcc = map(data, bootstrap_mcc, R = 100,
              ground = "ground_hemi", predicted = "pred_hemi",
              levels = c("LH", "RH"))

  ) %>%
  unnest(mcc)

summary_hemis <- results_hemis %>%
  group_by(method, hemis, hands, input) %>%
  summarize(
    n2 = n(),
    mcc_hemis_mean = weighted.mean(x = mcc_mean, w = n),
    mcc_hemis_sd = sqrt(Hmisc::wtd.var(mcc_mean, n)),
  ) %>%
  mutate(
    mcc_hemis_se = mcc_hemis_sd / sqrt(n2)
  )

results_hemis %>%
  group_by(hands) %>%
  summarize(
    min = min(mcc_mean),
    max = max(mcc_mean)
  )

### Plot results =====

png("plots/twoway.png", width = 6.5, height = 5, units = "in", res = 300)

ggplot(results_hemis, aes(x = input)) +
  geom_pointrange(aes(y = mcc_mean, ymin = mcc_ci_lo, ymax = mcc_ci_hi,
                      color = group),
                  position = position_dodge(width = 0.5),
                  size = 0.25, alpha = 0.75) +
  geom_point(data = summary_hemis,
              aes(x = input, y = mcc_hemis_mean)) +
  scale_x_discrete(limits = c("base", "base92", "osampl", "langtask"),
                   labels = c("HC", "HC92", "HCO", "LgHC")) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
  facet_grid(cols = vars(hands), rows = vars(method)) +
  theme_bw() +
  labs(x = "Method", y = "MCC", title = "Predict two-way hemisphere (only)",
       shape = "Input", linetype = "Input",
       color = "Test group") +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()


hemis_aov <- aov(mcc_mean ~ hands + method + input + group,
                 data = filter(results_hemis, !str_detect(input, "92")))

TukeyHSD(hemis_aov, c("hands", "method", "input"))

subset_check <- results_hemis %>%
  mutate(
    input2 = str_remove(input, "92"),
    subset = replace_na(subset, "none")
  ) %>%
  ungroup() %>%
  select(-data, -starts_with("mcc_ci"), -n) %>%
  select(method, group, hands, input, mcc_mean) %>%
  filter(
    hands != "all",
    input %in% c("base", "base92")
  ) %>%
  pivot_wider(names_from = c(hands, input), values_from = mcc_mean)

d00 <- t.test(subset_check$righty_base, subset_check$righty_base92, paired = TRUE)
d30 <- t.test(subset_check$righty2_base, subset_check$righty2_base92, paired = TRUE)

sin00 <- t.test(subset_check$lefty_base, subset_check$righty_base92, paired = TRUE)
sin30 <- t.test(subset_check$nonrighty_base, subset_check$righty2_base92, paired = TRUE)

t_tests <- tibble(test = list(d00, d30, sin00, sin30)) %>%
  mutate(
    effect = map_dbl(test, ~ .x$estimate[2] - .x$estimate[1]),
    dof = map_dbl(test, ~.x$parameter),
    t = map_dbl(test, ~round(.x$statistic, 2)),
    p = map_dbl(test, ~.x$p.value, 3)
  )

t_tests$p.adj <- p.adjust(t_tests$p, method = "fdr")

# subset_aov <- results_hemis %>%
#   mutate(
#     input2 = str_remove(input, "200"),
#     subset = replace_na(subset, "none")
#   ) %>%
#   filter(
#     !(hands %in% c("lefty", "nonrighty"))
#   ) %>%
#   aov(mcc_mean ~ subset + hands + method + input + group, data = .)
#
# TukeyHSD(subset_aov, c("hands", "method", "subset"))

## Predicting four-way handedness ====

results_hands1 <- results2 %>%
  filter(
    hemis == "all",
    outcome == "hand"
  ) %>%
  mutate(

    # Specify only evaluating hemi
    mcc = map(data, bootstrap_mcc, R = 100,
              ground = "ground_hand", predicted = "pred_hand",
              levels = c("lefty", "righty"),
              .progress = TRUE)

  ) %>%
  unnest(mcc)

results_hands2 <- results2 %>%
  filter(
    hemis == "all",
    outcome == "hand2"
  ) %>%
  mutate(

    # Specify only evaluating hemi
    mcc = map(data, bootstrap_mcc, R = 100,
              ground = "ground_hand", predicted = "pred_hand",
              levels = c("nonrighty", "righty2"),
              .progress = TRUE)

  ) %>%
  unnest(mcc)

results_hands <- bind_rows(results_hands1, results_hands2)

summary_hands <- results_hands %>%
  group_by(method, hemis, hands, input, outcome) %>%
  summarize(
    n2 = n(),
    mcc_hemis_mean = weighted.mean(x = mcc_mean, w = n),
    mcc_hemis_sd = sqrt(Hmisc::wtd.var(mcc_mean, n)),
  ) %>%
  mutate(
    mcc_hemis_se = mcc_hemis_sd / sqrt(n2)
  )

summary_hands %>%
  group_by(input, outcome) %>%
  summarize(
    min = min(mcc_hemis_mean),
    max = max(mcc_hemis_mean)
  ) %>%
  filter(
    outcome == "hand"
  )


# Poster plot ====

poster_data <- results_hemis %>%
  filter(
    hands %in% c("righty", "lefty")
  ) %>%
  mutate(
    method = toupper(method)
  )

summary_hemis_poster <- results_hemis %>%
  group_by(method, hemis, hands, input) %>%
  summarize(
    n2 = n(),
    mcc_hemis_mean = weighted.mean(x = mcc_mean, w = n),
    mcc_hemis_sd = sqrt(Hmisc::wtd.var(mcc_mean, n)),
  ) %>%
  mutate(
    mcc_hemis_se = mcc_hemis_sd / sqrt(n2)
  ) %>%
  filter(
    hands %in% c("righty", "lefty")
  ) %>%
  mutate(
    method = toupper(method)
  )


poster_results <- ggplot(NULL) +
  geom_point(data = poster_data, aes(x = input, y = mcc_mean, color = method),
             position = position_dodge(width = 0.5), alpha = 0.75,
             size = 5) +
  geom_pointrange(data = summary_hemis_poster,
                   shape = 21, size = 1.5,
                   position = position_dodge(width = 0.5),
                   aes(x = input, y = mcc_hemis_mean,
                       ymin = mcc_hemis_mean - mcc_hemis_se,
                       ymax = mcc_hemis_mean + mcc_hemis_se,
                       fill = method)) +
  scale_x_discrete(labels = c("Base", "n=92", "Lang.", "Over.")) +
  scale_y_continuous(limits = c(0.5, 1)) +
  facet_wrap(vars(hands)) +
  theme_bw() +
  theme(text = element_text(size = 24)) +
  labs(x = "Input", y = "MCC", color = "Method", fill = "Method")

ggsave("poster_plot.png", poster_results, height = 3.5, width = 9, units = "in")

summary_hands_poster <- results_hands1 %>%
  group_by(method, input) %>%
  summarize(
    n2 = n(),
    mcc_hemis_mean = weighted.mean(x = mcc_mean, w = n),
    mcc_hemis_sd = sqrt(Hmisc::wtd.var(mcc_mean, n)),
  ) %>%
  mutate(
    mcc_hemis_se = mcc_hemis_sd / sqrt(n2)
  ) %>%
  mutate(
    method = toupper(method)
  )

poster_plot_hands <- ggplot(NULL) +
  geom_point(data = mutate(results_hands1, method = toupper(method)),
             aes(x = input, y = mcc_mean, color = method),
             position = position_dodge(width = 0.5), alpha = 0.75,
             size = 5) +
  geom_pointrange(data = summary_hands_poster,
                  shape = 21, size = 1.5,
                  position = position_dodge(width = 0.5),
                  aes(x = input, y = mcc_hemis_mean,
                      ymin = mcc_hemis_mean - mcc_hemis_se,
                      ymax = mcc_hemis_mean + mcc_hemis_se,
                      fill = method)) +
  scale_x_discrete(labels = c("Base", "Full", "Lang.", "Over.",
                              "Trans.")) +
  scale_y_continuous(limits = c(NA, 1)) +
  theme_bw() +
  theme(text = element_text(size = 24)) +
  labs(x = "Input", y = "MCC", color = "Method", fill = "Method")

ggsave("poster_plot_hands.png", poster_plot_hands, height = 3.5, width = 7,
       units = "in")

png("plots/hands.png", width = 6.5, height = 5, units = "in", res = 300)

ggplot(results_hands, aes(x = input)) +
  geom_pointrange(aes(y = mcc_mean, ymin = mcc_ci_lo, ymax = mcc_ci_hi,
                      color = group),
                  position = position_dodge(width = 0.5),
                  size = 0.25, alpha = 0.75) +
  geom_hline(yintercept = c(0, 0.4), color = "red") +
  geom_point(data = summary_hands, aes(x = input, y = mcc_hemis_mean)) +
  scale_x_discrete(limits = c("base", "osampl", "langtask", "trans", "full"),
                    labels = c("HC", "HCO", "LgHC", "TC", "FC")) +
  facet_grid(cols = vars(outcome), rows = vars(method), scales = "free_x") +
  theme_bw() +
  labs(x = "Method", y = "MCC", title = "Predicting handedness",
       shape = "Input", linetype = "Input",
       color = "Test group") +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()

hands_aov <- aov(mcc_mean ~ method + input + group, data = results_hands)

TukeyHSD(hands_aov, which = c("input", "group"))

### Four-way handedness ====

results_4way <- results_big %>%
  filter(
    outcome %in% c("hand", "hand2"),
    input %in% c("base", "osampl", "langtask")
  ) %>%
  nest() %>%
  select(outcome, method, group, hands, hemis, input, data, n) %>%
  mutate(

    lefty_data = map(data, filter, str_detect(ground, "lefty|nonrighty")),
    righty_data = map(data, filter, str_detect(ground, "^righty")),

  )

results_4way_hand1 <- results_4way %>%
  filter(
    outcome == "hand"
  ) %>%
  mutate(

    all = map(data, ~bootstrap_mcc(.x, R = 100,
                                   levels = c("leftyLH", "leftyRH",
                                              "rightyLH", "rightyRH"))),

    lefty = map(lefty_data, ~bootstrap_mcc(.x, R = 100,
                                     levels = c("leftyLH", "leftyRH",
                                                "rightyLH", "rightyRH"))),

    righty = map(righty_data, ~bootstrap_mcc(.x, R = 100,
                                           levels = c("leftyLH", "leftyRH",
                                                      "rightyLH", "rightyRH"))),

  ) %>%
  unnest(c(all, lefty, righty), names_sep = "_")

results_4way_hand2 <- results_4way %>%
  filter(
    outcome == "hand2"
  ) %>%
  mutate(

    all = map(data, ~bootstrap_mcc(.x, R = 100,
                                   levels = c("nonrightyLH", "nonrightyRH",
                                              "righty2LH", "righty2RH"))),

    nonrighty = map(lefty_data,
                ~bootstrap_mcc(.x, R = 100,
                               levels = c("nonrightyLH", "nonrightyRH",
                                          "righty2LH", "righty2RH"))),

    righty2 = map(righty_data,
                 ~bootstrap_mcc(.x, R = 100,
                                levels = c("nonrightyLH", "nonrightyRH",
                                           "righty2LH", "righty2RH"))),
  ) %>%
  unnest(c(all, nonrighty, righty2), names_sep = "_")

results_4way_2 <- bind_rows(results_4way_hand1, results_4way_hand2) %>%
  select(outcome, method, group, input, contains("mcc")) %>%
  pivot_longer(contains("mcc")) %>%
  separate_wider_delim(name, delim = "_", names = c("hgroup", NA, "name"),
                       too_many = "merge") %>%
  pivot_wider(names_from = name, values_from = value) %>%
  filter(
    !is.na(mean)
  )

range(results_4way_2$mean[results_4way_2$hgroup == "all"])

ggplot(results_4way_2, aes(x = input)) +
  geom_pointrange(aes(y = mean, ymin = ci_lo, ymax = ci_hi,
                      color = group),
                  position = position_dodge(width = 0.5),
                  size = 0.25, alpha = 0.5) +
  geom_hline(yintercept = 0, color = "red") +
  scale_x_discrete(limits = c("base", "osampl", "langtask"),
                     labels = c("HC", "HCO", "LgHC")) +
  scale_y_continuous(limits = c(NA, NA)) +
  facet_grid(rows = vars(method), cols = vars(hgroup)) +
  theme_bw() +
  labs(x = "Method", y = "MCC", color = "Handedness group",
       title = "Predicting four-way distinction") +
  theme(legend.position = "bottom") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")


fourway_aov <- aov(mean ~ method + input + hgroup + group,
                   data = results_4way_2)

TukeyHSD(fourway_aov, which = c("input", "hgroup", "group"))
