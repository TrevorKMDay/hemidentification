library(tidyverse)

select <- dplyr::select

setwd("~/Google Drive/My Drive/Projects/hemisphere_fingerprinting/code/modeling/")

# demographic data ====

hands <- read_csv("../../data/restricted_data.csv", show_col_types = FALSE) %>%
  rename(
    EHI = Handedness
  ) %>%
  mutate(
    sub = paste0("sub-", Subject),
    handedness = if_else(EHI > 0, "righty", "lefty"),
    handedness2 = if_else(EHI >= 30, "righty2", "nonrighty")
  ) %>%
  select(sub, EHI, starts_with("handedness"))

ggplot(hands, aes(x = EHI)) +
  geom_histogram(binwidth = 10, boundary = 5) +
  scale_x_continuous(breaks = seq(-100, 100, 20)) +
  theme_bw()

ggplot(hands, aes(x = EHI)) +
  geom_density() +
  scale_x_continuous(breaks = seq(-100, 100, 20)) +
  theme_bw()


hands_var <- function(ehi_vec, thresh) {

  left <- ehi_vec[ehi_vec <= thresh]
  right <- ehi_vec[ehi_vec > thresh]

  return(tibble(left = sd(left), right = sd(right)))

}

hands_var_by_thresh <- tibble(thresh = seq(-90, 90, by = 10)) %>%
  mutate(
    sd = map(thresh, ~hands_var(hands$EHI, .x))
  ) %>%
  unnest(sd) %>%
  pivot_longer(-thresh)

ggplot(hands_var_by_thresh, aes(x = thresh, y = value)) +
  geom_line(aes(color = name))

table(hands$EHI < -75)
