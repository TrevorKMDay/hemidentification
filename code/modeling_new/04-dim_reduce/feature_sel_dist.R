library(tidyverse)
library(segmented)

select <- dplyr::select

home <- "/Users/tkmd/MyDrive/Projects/hemisphere_fingerprinting/code/modeling_new/04-dim_reduce"
setwd(home)

all_in_files <- list.files("all_in/", pattern = "results.csv", full.names = TRUE)

results0 <- tibble(f = all_in_files) %>%
  mutate(
    data = map(f, read_csv, show_col_types = FALSE)
  )

results <- results0 %>%
  mutate(
    bn = basename(f)
  ) %>%
  separate_wider_delim(bn, delim = "_", names = c("input", "k", NA, NA, NA)) %>%
  unnest(data)  %>%
  select(sub, EHI, input, hemi, k, LD1) %>%
  pivot_wider(names_from = hemi, values_from = LD1) %>%
  group_by(input) %>%
  mutate( 

    dist = RH - LH,
    distZ = scale(dist, center = FALSE)[, 1],
    
    k = str_remove(k, "k-") %>%
      as.numeric(),
    k_log = log(k, 2)

  )

lefties <- unique(results$sub[results$input == "data-hclhd"])

best_model <- results %>%
  filter(
    k == 16110,
    input %in% c("data-hclhd", "data-hcrhd", "data-hemiconnectome")
  ) %>%
  mutate(
    hand = if_else(sub %in% lefties, "lefty", "righty"),
    input = replace_values(
      input,
      "data-hclhd" ~ "Sinistrals only",
      "data-hcrhd" ~ "Dextrals only",
      "data-hemiconnectome" ~ "All participants"
    ) 
  )

best_model_cors <- best_model %>%
  group_by(input) %>%
  nest() %>%
  mutate(
    r = map(data, ~Hmisc::rcorr(.x$LH, .x$RH))
  )

ggplot(best_model, aes(x = LH, y = RH)) +
  geom_point() +
  facet_wrap(vars(input)) +
  theme_bw()

ggplot(best_model, aes(x = EHI, y = dist)) +
  geom_point(alpha = 0.25) +
  geom_smooth(method = "lm") +
  geom_vline(xintercept = 0, color = "red") +
  facet_wrap(vars(input), scales = "free") +
  theme_bw() +
  labs(x = "EHI", y = "Hemisphere distance")

ggsave("../plots/ehi_distance.png", width = 6.5, height = 3.5)

lm(dist ~ EHI, data = best_model_cors$data[[1]]) %>%
  summary()

lm(dist ~ EHI, data = best_model_cors$data[[2]]) %>%
  summary()

lm(dist ~ EHI, data = best_model_cors$data[[3]]) %>%
  summary()

# Show how distance changes with number of features ====

ggplot(results, aes(x = k_log, y = dist, color = input)) +
  geom_point(alpha = 0.1) +
  geom_line(aes(group = interaction(sub, input)), alpha = 0.1) +
  scale_x_continuous(breaks = 0:14, labels = 2^(0:14)) +
  theme_bw()
 
ggplot(filter(results, k == 16110), aes(x = EHI, y = dist)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(vars(input)) +
  theme_bw()
 
ggplot(filter(results, input == "data-hemiconnectome"), aes(x = EHI, y = dist)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(vars(k), scales = "free") +
  theme_bw()

variance <- results %>%
  group_by(input, k, k_log) %>%
  summarize(
    var = var(dist),
    sd = sd(dist),
    mean = mean(dist), 
  ) %>%
  mutate(
    coef_var = sd /mean
  )

ggplot(variance, aes(x = k_log, y = var, color = input)) +
  geom_point() +
  geom_line() +
  theme_bw()

ggplot(variance, aes(x = k_log, y = sd, color = input)) +
  geom_point() +
  geom_line() +
  theme_bw()

ggplot(variance, aes(x = k_log, y = coef_var, color = input)) +
  geom_point() +
  geom_line() +
  theme_bw()

# Segmented regression ====

best_model_all <- best_model %>%
  filter(
    input == "All participants"
  )

lm0 <- lm(dist ~ 1, data = best_model_all)
lm1 <- lm(dist ~ EHI, data = best_model_all)

lm1_seg1 <- segmented(lm1, seg.Z = ~EHI, npsi = 1)
anova(lm0, lm1, lm1_seg1) # OK
AIC(lm0, lm1, lm1_seg1) # OK

lm1_seg2 <- segmented(lm1, seg.Z = ~EHI, npsi = 2)
anova(lm0, lm1, lm1_seg1, lm1_seg2) # OK
AIC(lm0, lm1, lm1_seg1, lm1_seg2) # OK

lm1_seg3 <- segmented(lm1, seg.Z = ~EHI, npsi = 3)
anova(lm0, lm1, lm1_seg1, lm1_seg2, lm1_seg3) # OK
AIC(lm0, lm1, lm1_seg1, lm1_seg2, lm1_seg3) # OK

best_model_segmented <- best_model_all %>%
  mutate(
    seg1_groups = if_else(EHI > lm1_seg1$psi[2], "g1", "g2")
  )

ggplot(best_model_segmented, aes(x = EHI, y = dist)) +
  geom_jitter(width = 0.25, height = 0, alpha = 0.5)  + 
  geom_smooth(aes(group = seg1_groups), method = "lm") +
  geom_vline(xintercept = 0.5, color = "red") +
  theme_bw() +
  labs(y = "Hemisphere distance")

ggsave("plots/ehi_distance_1bp.png", width = 6.5, height = 3.5)
