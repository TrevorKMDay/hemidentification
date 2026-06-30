library(tidyverse)
library(qs2)
library(effectsize)

setwd("~/MyDrive/Projects/hemisphere_fingerprinting/code/modeling/")
hc <- qs_read("inputs/hemiconnectome.rds")

hc2 <- hc %>%
  select(sub, hemi, handedness, EHI, PF_PSL, IFSa_PSL, `6ma_PSL`)

hc3 <- hc2 %>%
  pivot_longer(ends_with("_PSL"), names_to = "cxn")

ggplot(hc3, aes(value)) +
  geom_density(aes(fill = hemi)) +
  facet_wrap(vars(cxn), scales = "free") +
  theme_bw()

ggplot(hc3, aes(x = EHI, y = value, color = hemi)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "lm") +
  geom_hline(yintercept = 0, color = "red") +
  facet_wrap(vars(cxn), scales = "free") +
  theme_bw()

hc4 <- hc3 %>%
  pivot_wider(names_from = hemi) %>%
  group_by(cxn) %>%
  nest() %>%
  mutate(
    d = map(data, ~cohens_d(x = .x$LH, y = .x$RH, paired = TRUE))
  )

hc3 %>%
  pivot_wider(names_from = hemi) %>%
  ggplot(aes(x = LH, y = RH)) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1) +
    facet_wrap(vars(cxn)) +
    coord_equal()
