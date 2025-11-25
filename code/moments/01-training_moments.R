library(tidyverse)
library(qs)
library(moments)
library(esc)
library(PearsonDS)

setwd("~/MyDrive/Projects/hemisphere_fingerprinting/code/moments/")

# Load data ====

# Fisher-Z data
hc0 <- qread("../modeling/inputs/hemiconnectome.rds")

# Z(FisherZ(r))
hcz0 <- qread("../modeling/inputs/hemiconnectome_Z.rds")

id_cols <- c("sub", "gender", "age_group", "hemi", "handedness")

# Regular data ====

hc <- hc0 %>%
  select(all_of(id_cols), contains("_")) %>%
  pivot_longer(-all_of(id_cols), values_to = "r")

hcZ <- hcz0 %>%
  select(all_of(id_cols), contains("_")) %>%
  pivot_longer(-all_of(id_cols), values_to = "rZ")

hc_moments <- hc %>%
  group_by(across(all_of(id_cols))) %>%
  summarize(
    mean = mean(r),
    var = var(r),
    skew = skewness(r),
    kurt = kurtosis(r)
  )

hcZ_moments <- hcZ %>%
  group_by(across(all_of(id_cols))) %>%
  summarize(
    mean = mean(rZ),
    var = var(rZ),
    skew = skewness(rZ),
    kurt = kurtosis(rZ)
  )

hc_moments_long <- hc_moments %>%
  pivot_longer(c(mean, var, skew, kurt), names_to = "moment")

hcZ_moments_long <- hcZ_moments %>%
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


x <- seq(-1, 2, by = 0.1)
x2 <- seq(-3.5, 3.5, by = 0.1)

mean_curves <- hc_moments %>%
  group_by(hemi) %>%
  summarize(
    across(c(mean, var, skew, kurt), mean)
  )

mean_Z_curve <- hcZ_moments %>%
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

hcZ_plot <- hcZ_moments %>%
  rowwise() %>%
  mutate(
    curve = list(tibble(x = x2,
                        y = dpearson(x2,
                                     moments = c(mean = mean,
                                                 variance = var,
                                                 skewness = skew,
                                                 kurtosis = kurt))))
  ) %>%
  ungroup() %>%
  select(sub, hemi, handedness, curve) %>%
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

mean_Z_plot <- mean_Z_curve %>%
  rowwise() %>%
  mutate(
    curve = list(tibble(x = x2,
                        y = dpearson(x2,
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
  theme_bw() +
  labs(title = "Fisher-Z transformed data", x = "FisherZ(r)", y = "Density")

ggplot(hcZ_plot, aes(x = x, y = y, color = hemi)) +
  geom_line(aes(group = sub), alpha = 0.1) +
  geom_line(
    data = mean_Z_plot,
    color = "black"
  ) +
  facet_wrap(vars(hemi)) +
  theme_bw() +
  labs(title = "Normalized, Fisher-Z transformed data", x = "FisherZ(r)",
       y = "Density")

# ggplot(mean_plot, aes(x = x, y = y, color = hemi)) +
#   geom_line()

# Lefties ===

# FisherZ(r)
lefties0 <- read_csv("~/Projects/apply_hemiconnectome/lc1_stat-rfZ.csv",
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

# Z(FisherZ(r))
leftiesZ0 <- read_csv("~/Projects/apply_hemiconnectome/lc1_stat-rfZZ.csv",
                      show_col_types = FALSE)

leftiesZ <- leftiesZ0 %>%
  mutate(
    bn = basename(f),
    sub = str_extract(bn, "sub-LC...")
  ) %>%
  filter(
    !str_detect(bn, "run-1")
  ) %>%
  select(sub, hemi, everything(), -f)

leftiesZ_moments <- leftiesZ %>%
  pivot_longer(-c(bn, sub, hemi), values_to = "r") %>%
  group_by(sub, hemi) %>%
  summarize(
    mean = mean(r),
    var = var(r),
    skew = skewness(r),
    kurt = kurtosis(r)
  )

leftiesZ_moments_long <- leftiesZ_moments %>%
  pivot_longer(c(mean, var, skew, kurt), names_to = "moment")


leftiesZ_curves <- leftiesZ_moments %>%
  rowwise() %>%
  mutate(
    curve = list(tibble(x = x2,
                        y = dpearson(x2,
                                     moments = c(mean = mean,
                                                 variance = var,
                                                 skewness = skew,
                                                 kurtosis = kurt))))
  ) %>%
  ungroup() %>%
  select(sub, hemi, curve) %>%
  unnest(curve)

ggplot(leftiesZ_curves, aes(x = x, y = y)) +
  geom_line(aes(color = sub)) +
  facet_wrap(vars(hemi)) +
  theme_bw() +
  labs(x = "Z(FisherZ(r)) for lefties")


# Combined plot ====

ggplot(hc_plot, aes(x = x, y = y)) +
  geom_line(aes(group = sub), alpha = 0.0075) +
  geom_line(data = lefties_curves, aes(color = sub)) +
  facet_wrap(vars(hemi)) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(title = "Lefty athletes overlaid on HCP", x = "FisherZ(r)",
       y = "Density")

ggplot(hcZ_plot, aes(x = x, y = y)) +
  geom_line(aes(group = sub), alpha = 0.0075) +
  geom_line(data = leftiesZ_curves, aes(color = sub)) +
  facet_wrap(vars(hemi)) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(title = "Lefty athletes overlaid on HCP, normalized",
       x = "Z(FisherZ(r))", y = "Density")

z_hcprighties <- ggplot(filter(hcZ_plot, handedness == "righty"),
                        aes(x = x, y = y)) +
  geom_line(aes(group = sub), alpha = 0.0075) +
  geom_line(data = leftiesZ_curves, aes(color = sub)) +
  facet_wrap(vars(hemi)) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(title = "Lefty athletes overlaid on HCP righties, normalized",
       x = "Z(FisherZ(r))", y = "Density")

z_hcplefties <- ggplot(filter(hcZ_plot, handedness == "lefty"),
                       aes(x = x, y = y)) +
  geom_line(aes(group = sub), alpha = 0.1) +
  geom_line(data = leftiesZ_curves, aes(color = sub)) +
  facet_wrap(vars(hemi)) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(title = "Lefty athletes overlaid on HCP lefties, normalized",
       x = "Z(FisherZ(r))", y = "Density")

library(patchwork)

z_hcprighties / z_hcplefties +
  plot_layout(axes = "collect", guides = "collect") &
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

# Short HC ====

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
    aes(group = sub),
    alpha = 0.
  ) +
  geom_line(
    data = lefties_curves,
    aes(color = sub),
  ) +
  facet_wrap(vars(hemi)) +
  theme_bw() +
  theme(legend.position = "none")

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

# Compare Z(FisherZ) ====

test <- tibble(r = hc$r[hc$sub == "sub-100206"],
               rZ = hcZ$rZ[hcZ$sub == "sub-100206"])

yhat <- filter(hc_plot, sub == "sub-100206")
yhat2 <- filter(hcZ_plot, sub == "sub-100206")

ggplot(NULL) +
  geom_density(data = test, aes(r), color = "red") +
  geom_density(data = test, aes(rZ), color = "orange") +
  geom_line(
    data = yhat2,
    aes(x = x, y = y, color = hemi)
  ) +
  theme_bw()

ggplot(hcZ_plot, aes(x = x, y = y)) +
  geom_line(aes(group = sub), alpha = 0.0075) +
  geom_line(data = leftiesZ_curves, aes(color = sub)) +
  facet_wrap(vars(hemi)) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(title = "Lefty athletes overlaid on HCP", x = "FisherZ(r)",
       y = "Density")

ggplot(hcZ_moments_long, aes(x = hemi, fill = hemi, y = value)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_point(
    data = leftiesZ_moments_long,
    aes(shape = sub),
    fill = "red", size = 2.5
  ) +
  scale_shape_manual(values = c(21, 22, 23, 24)) +
  facet_wrap(vars(moment), scales = "free") +
  theme_bw()


# Short HC ====

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
    aes(group = sub),
    alpha = 0.
  ) +
  geom_line(
    data = lefties_curves,
    aes(color = sub),
  ) +
  facet_wrap(vars(hemi)) +
  theme_bw() +
  theme(legend.position = "none")

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
