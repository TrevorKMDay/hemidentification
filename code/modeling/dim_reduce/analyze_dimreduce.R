library(tidyverse)

setwd("~/MyDrive/Projects/hemisphere_fingerprinting/code/modeling/dim_reduce/")

# Scores ====

score_files <- list.files(pattern = "*score*")

scores0 <- tibble(f = score_files) %>%
  mutate(
    data = map(f, read_csv, show_col_types = FALSE)
  )

scores <- scores0 %>%
  unnest(data) %>%
  mutate(
    group = str_extract(f, "group-.") %>%
      str_remove("group-")
  ) %>%
  select(group, feature, score)

scores2 <- scores %>%
  group_by(feature) %>%
  summarize(
    n = n(),
    mean_score = mean(score),
    sd_score = sd(score)
  ) %>%
  arrange(
    desc(mean_score)
  ) %>%
  mutate(
    se_score = sd_score / sqrt(n),
    rank = row_number(),
    label = if_else(rank <= 10, feature, NA_character_)
  )

ggplot(scores2, aes(x = rank, y = mean_score)) +
  geom_pointrange(aes(ymin = mean_score - se_score,
                      ymax = mean_score + se_score)) +
  theme_bw()

str_subset(scores2$feature, "_", negate = TRUE)

scores3 <- scores2 %>%
  filter(
    rank <= 25
  ) %>%
  select(feature, mean_score, se_score)

clipr::write_clip(scores3)

top512 <- scores2$feature[scores2$rank <= 512]

# accuracy

result_files <- list.files(pattern = ".*result.*")

results0 <- tibble(f = result_files) %>%
  mutate(
    data = map(f, read_csv, show_col_types = FALSE)
  )

results <- results0 %>%
  mutate(

    group = str_extract(f, "group-.") %>%
      str_remove("group-"),

    k = str_extract(f, "k-[0-9]*") %>%
      str_remove("k-") %>%
      as.numeric(),

  ) %>%
  select(group, k, data) %>%
  mutate(
    n = map_int(data, nrow),
    correct = map_int(data, ~sum(.x$hemi == .x$predicted)),
    p = correct / n
  )

ggplot(results, aes(x = k, y = p, color = group)) +
  geom_point() +
  geom_line() +
  scale_y_continuous(limits = c(0.85, NA)) +
  theme_bw()

ggplot(results, aes(x = log(k), y = p, color = group)) +
  geom_point() +
  geom_line() +
  scale_y_continuous(limits = c(0.85, NA)) +
  theme_bw()

results_summary <- results %>%
  group_by(k) %>%
  summarize(
    n = n(),
    mean = mean(p),
    sd_p = sd(p)
  ) %>%
  mutate(
    se_p = sd_p / sqrt(n)
  )

ggplot(results_summary, aes(x = k, y = mean)) +
  geom_pointrange(aes(ymin = mean - se_p, ymax = mean + se_p)) +
  geom_line() +
  scale_y_continuous(limits = c(NA, NA)) +
  theme_bw()

# Which hemi is >

library(qs2)
library(effectsize)

hc <- qs_read("../inputs/hemiconnectome.rds")

hc2 <- hc %>%
  select(-gender, -age_group, -starts_with('handedness'), -starts_with("class"),
         -EHI, -group)

hc3 <- hc2 %>%
  pivot_longer(-c(sub, hemi), names_to = "cxn")

# ggplot(hc3, aes(value)) +
#   geom_density(aes(fill = hemi)) +
#   facet_wrap(vars(cxn), scales = "free") +
#   theme_bw()

# ggplot(hc3, aes(x = EHI, y = value, color = hemi)) +
#   geom_point(alpha = 0.1) +
#   geom_smooth(method = "lm") +
#   geom_hline(yintercept = 0, color = "red") +
#   facet_wrap(vars(cxn), scales = "free") +
#   theme_bw()

hc4 <- hc3 %>%
  pivot_wider(names_from = hemi) %>%
  group_by(cxn) %>%
  nest() %>%
  mutate(
    d = map(data,
            ~cohens_d(x = .x$LH, y = .x$RH, paired = TRUE, verbose = FALSE),
            .progress = TRUE)
  )

hc5 <- hc4 %>%
  mutate(
    effsize = map_dbl(d, ~ .x$Cohens_d)
  ) %>%
  select(cxn, effsize)

rank_diff <- left_join(scores2, hc5, by = join_by(feature == cxn))

ggplot(rank_diff, aes(x = effsize, y = mean_score, color = rank > 512)) +
  geom_point(alpha = 0.1) +
  labs(x = "RH>LH --- LH>RH") +
  theme_bw()
