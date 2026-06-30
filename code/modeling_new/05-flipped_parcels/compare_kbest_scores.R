library(tidyverse)

mirrored_score_files <- list.files(pattern = "selKbestscores")

scores <- tibble(f = mirrored_score_files) %>%
  mutate(
    scores = map(f, read_csv, show_col_types = FALSE),
    model = basename(f) %>%
      str_remove_all("_selKbestscores.csv|data-")
  ) %>% 
  unnest(scores) %>%
  select(-f) %>%
  pivot_wider(names_from = model, values_from = score)

score_s_corr <- corrr::correlate(scores, method = "s")
clipr::write_clip(score_s_corr)

lefty_scores <- read_csv("../lefty_scores.csv")
right_scores <- read_csv("../righty_scores.csv")

all_scores <- scores %>%
  pivot_longer(-feature) %>%
  separate_wider_delim(name, ".", names = c(NA, "flip", "handedness")) %>%
  pivot_wider(names_from = handedness) %>%
  left_join(lefty_scores) %>%
  right_join(right_scores) %>%
  group_by(flip) %>%
  mutate(
    across(c(LHD, RHD, ends_with("score")), ~rank(-.x), .names = "{col}_rank")
  )


ggplot(all_scores, aes(x = scale(righty_score), y = scale(RHD), color = flip)) +
  geom_point(alpha = 0.25) +
  geom_smooth(method = "lm") +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  theme_bw() +
  labs(x = "Regular score", y = "Mirrored score", title = "Righties")

ggplot(all_scores, aes(x = scale(lefty_score), y = scale(LHD), color = flip)) +
  geom_point(alpha = 0.25) +
  geom_smooth(method = "lm") +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  theme_bw() +
  labs(x = "Regular score", y = "Mirrored score", title = "Lefties")

ggplot(all_scores, aes(x = righty_score_rank, y = RHD_rank, color = flip)) +
  geom_point(alpha = 0.25) +
  geom_smooth(method = "lm") +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  theme_bw()

ggplot(all_scores, aes(x = lefty_score_rank, y = LHD_rank, color = flip)) +
  geom_point(alpha = 0.25) +
  geom_smooth(method = "lm") +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  theme_bw()

all_scores_corr <- all_scores %>%
  group_by(flip) %>%
  summarize(
    righty_r = cor(RHD, righty_score, method = "s"),
    lefty_r = cor(LHD, lefty_score, method = "s")
  )

all_scores_top10pct <- all_scores %>%
  filter(
    LHD_rank < 1611 | RHD_rank < 1611 | lefty_score_rank < 1611 | 
      righty_score_rank < 1611
  )

ggplot(all_scores, aes(x = righty_score, y =RHD, color = flip)) +
  geom_point(alpha = 0.25) +
  geom_smooth(method = "lm") +
  facet_wrap(vars(flip)) +
  theme_bw()

ggplot(all_scores, aes(x = lefty_score, y = LHD, color = flip)) +
  geom_point(alpha = 0.25) +
  geom_smooth(method = "lm") +
  facet_wrap(vars(flip)) +
  theme_bw()


all_scores_top10pct_corr <- all_scores_top10pct %>%
  group_by(flip) %>%
  summarize(
    righty_r = cor(RHD, righty_score, method = "p"),
    lefty_r = cor(LHD, lefty_score, method = "p")
  )

# Difference and compare to asymmetry

asym_index <- read_csv("~/MyDrive/Projects/parcel_overlap/glasser_LR_overlap.csv") %>%
  select(label, dice)

all_scores2 <- all_scores %>%
  select(-ends_with("rank")) %>%
  separate_wider_delim(feature, delim = "_", names = c(NA, "ROI1", "ROI2")) %>%
  left_join(
    asym_index, by = join_by(ROI1 == label)
  ) %>%
  left_join(
    asym_index, by = join_by(ROI2 == label), suffix = c("1", "2")
  ) %>% 
  mutate(
    across(c(LHD, RHD, lefty_score, righty_score), ~scale(.x)[,1]),
    lefty_diff = lefty_score - LHD,
    righty_diff = righty_score - RHD,
    avg_dice = (dice1 + dice2) / 2
  ) 

all_scores_plotme <- all_scores2 %>%
  select(ROI1, ROI2, flip, avg_dice, ends_with("diff")) %>%
  pivot_longer(ends_with("diff"), values_to = "diff")

ggplot(all_scores_plotme, aes(x = avg_dice, y = diff)) + 
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm") +
  facet_grid(rows = vars(flip), cols = vars(name)) +
  theme_bw()

top_difference <- all_scores2 %>%
  filter(
    abs(lefty_diff) > 9 | abs(righty_diff) > 9
  ) %>%
  select(ROI1, ROI2, avg_dice, ends_with("_diff"))

clipr::write_clip(top_difference)
