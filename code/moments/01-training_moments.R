library(tidyverse)
library(qs)
library(moments)

setwd("~/MyDrive/Projects/hemisphere_fingerprinting/code/moments/")

hc0 <- qread("../modeling/inputs/hemiconnectome.rds")

id_cols <- c("sub", "gender", "age_group", "hemi", "handedness")

hc <- hc0 %>%
  select(all_of(id_cols), contains("_")) %>%
  pivot_longer(-all_of(id_cols), values_to = "r")

hc_moments <- hc %>%
  group_by(across(all_of(id_cols))) %>%
  summarize(
    mean = mean(r),
    var = var(r),
    skew = skewness(r),
    kurt = kurtosis(r)
  )

hc_moments_long <- hc_moments %>%
  pivot_longer(c(mean, var, skew, kurt), names_to = "moment")


hc_moments_hemi <- hc_moments_long %>%
  group_by(moment) %>%
  nest() %>%
  mutate(
    test = map(data, ~lm(value ~ hemi, data = .x)),

    coef = map_dbl(test, ~coef(.x)[2]),
    sdY = map_dbl(data, ~sd(.x$value)),
    n = map_int(data, nrow),

    eff_size = map(test, ~esc_B(coef, sdY, grp1n = n, grp2n = n)),
    d = map_dbl(eff_size, ~ .x$es)

  )

# There is an effect of hemisphere:

# Mean: RH > LH
# Var:  RH > LH
# Skew: RH < LH
# Kurt: RH < LH

hc_moments_hand <- hc_moments_long %>%
  group_by(moment) %>%
  nest() %>%
  mutate(
    test = map(data, ~lm(value ~ handedness + hemi, data = .x)),

    # coef = map_dbl(test, ~coef(.x)[2]),
    # sdY = map_dbl(data, ~sd(.x$value)),
    # n = map_int(data, nrow),
    #
    # eff_size = map(test, ~esc_B(coef, sdY, grp1n = n, grp2n = n)),
    # d = map_dbl(eff_size, ~ .x$es)

  )

# There is no effect of handedness

# plot me

library(PearsonDS)

x <- seq(-1, 2, by = 0.05)

mean_curves <- hc_moments %>%
  group_by(hemi) %>%
  summarize(
    across(c(mean, var, skew, kurt), mean)
  )

hc_plot <- hc_moments %>%
  rowwise() %>%
  mutate(
    curve = list(tibble(x = x,
                        y = dpearson(x,
                                     moments = c(mean = mean,
                                                 variance = var,
                                                 skewness = skew,
                                                 kurtosis = kurt))))
  ) %>%
  ungroup() %>%
  select(sub, hemi, curve) %>%
  unnest(curve)

mean_plot <- mean_curves %>%
  rowwise() %>%
  mutate(
    curve = list(tibble(x = x,
                        y = dpearson(x,
                                     moments = c(mean = mean,
                                                 variance = var,
                                                 skewness = skew,
                                                 kurtosis = kurt))))
  ) %>%
  unnest(curve)

ggplot(hc_plot, aes(x = x, y = y, color = hemi)) +
  geom_line(aes(group = sub), alpha = 0.1) +
  geom_line(
    data = mean_plot,
    color = "black"
  ) +
  facet_wrap(vars(hemi)) +
  theme_bw()

ggplot(mean_plot, aes(x = x, y = y, color = hemi)) +
  geom_line()

# Lefties ===

lefties0 <- read_csv("~/Projects/apply_hemiconnectome/lc1.csv",
                     show_col_types = FALSE)

lefties <- lefties0 %>%
  mutate(
    bn = basename(f),
    sub = str_extract(bn, "sub-LC...")
  ) %>%
  filter(
    !str_detect(bn, "run-1")
  ) %>%
  select(sub, hemi, everything(), -f)

lefties_moments <- lefties %>%
  pivot_longer(-c(bn, sub, hemi), values_to = "r") %>%
  group_by(sub, hemi) %>%
  summarize(
    mean = mean(r),
    var = var(r),
    skew = skewness(r),
    kurt = kurtosis(r)
  )

lefties_moments_long <- lefties_moments %>%
  pivot_longer(c(mean, var, skew, kurt), names_to = "moment")

lefties_curves <- lefties_moments %>%
  rowwise() %>%
  mutate(
    curve = list(tibble(x = x,
                        y = dpearson(x,
                                     moments = c(mean = mean,
                                                 variance = var,
                                                 skewness = skew,
                                                 kurtosis = kurt))))
  ) %>%
  ungroup() %>%
  select(sub, hemi, curve) %>%
  unnest(curve)

ggplot(lefties_curves, aes(x = x, y = y)) +
  geom_line(aes(color = sub)) +
  facet_wrap(vars(hemi)) +
  theme_bw()

# Combined plot ====

ggplot(hc_plot, aes(x = x, y = y)) +
  geom_line(aes(group = sub), alpha = 0.0075) +
  geom_line(data = lefties_curves, aes(color = sub)) +
  facet_wrap(vars(hemi)) +
  theme_bw() +
  theme(legend.position = "bottom")

ggplot(hc_moments_long, aes(x = hemi, fill = hemi, y = value)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_point(
    data = lefties_moments_long,
    aes(shape = sub),
    fill = "red", size = 2.5
  ) +
  scale_shape_manual(values = c(21, 22, 23, 24)) +
  facet_wrap(vars(moment), scales = "free") +
  theme_bw()

# Short hc ====

short_hc0 <- read_csv("../../data/shortened_pconns/short_hemiconnectome.csv",
                      show_col_types = FALSE)

short_hc <- short_hc0

shorthc_moments <- short_hc %>%
  pivot_longer(-c(sub, hemi), values_to = "r") %>%
  group_by(sub, hemi) %>%
  summarize(
    mean = mean(r),
    var = var(r),
    skew = skewness(r),
    kurt = kurtosis(r)
  )

shorthc_moments_long <- shorthc_moments %>%
  pivot_longer(c(mean, var, skew, kurt), names_to = "moment")

shorthc_curves <- shorthc_moments %>%
  rowwise() %>%
  mutate(
    curve = list(tibble(x = x,
                        y = dpearson(x,
                                     moments = c(mean = mean,
                                                 variance = var,
                                                 skewness = skew,
                                                 kurtosis = kurt))))
  ) %>%
  ungroup() %>%
  select(sub, hemi, curve) %>%
  unnest(curve)

ggplot(hc_plot, aes(x = x, y = y)) +
  geom_line(aes(group = sub), alpha = 0.0075) +
  geom_line(
    data = shorthc_curves,
    aes(color = sub),
    alpha = 0.5
  ) +
  geom_line(
    data = lefties_curves,
    aes(group = sub),
    color = "red"
    ) +
  facet_wrap(vars(hemi)) +
  theme_bw() +
  theme(legend.position = "none")

ggplot(hc_moments_long, aes(x = hemi, fill = hemi, y = value)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_point(
    data = shorthc_moments_long,
    fill = "red", size = 2.5
  ) +
  geom_point(
    data = lefties_moments_long,
    aes(shape = sub),
    fill = "red", size = 2.5
  ) +
  scale_shape_manual(values = c(21, 22, 23, 24)) +
  facet_wrap(vars(moment), scales = "free") +
  theme_bw()

ggplot(shorthc_moments_long, aes(x = hemi, fill = hemi, y = value)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_point(
    data = lefties_moments_long,
    aes(shape = sub),
    fill = "red", size = 2.5
  ) +
  scale_shape_manual(values = c(21, 22, 23, 24)) +
  facet_wrap(vars(moment), scales = "free") +
  theme_bw()
