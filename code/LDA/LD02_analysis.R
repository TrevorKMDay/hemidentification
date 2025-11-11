library(tidyverse)
library(igraph)

# Setup ====

select <- dplyr::select

setwd("~/MyDrive/Projects/hemisphere_fingerprinting/code/LDA")

# Load data ====

hc <- qs::qread("modeling/inputs/hemiconnectome.rds")

table(hc$handedness) / 2

files <- list.files("modeling/results/base", "*_xform.csv", full.names = TRUE)
p_files <- list.files("modeling/results/base", "*_predicted.csv",
                     full.names = TRUE)

pred0 <- tibble(f = p_files) %>%
  mutate(
    data = map(f, read_csv, show_col_types = FALSE)
  ) %>%
  unnest(data)

xform0 <- tibble(f = files) %>%
  mutate(
    data = map(f, read_csv, show_col_types = FALSE)
  )

xform1 <- xform0 %>%
  unnest(data) %>%
  select(-f) %>%
  rename(
    LD1 = `0`,
    LD2 = `1`,
    LD3 = `2`,
    ground = class
  )

# This produces strong separation in classes, is this from the self train/test
# sample?

# xform1 <- read_csv("modeling/lda_full/lda_class_xform.csv",
#                    show_col_types = FALSE) %>%
#   select(-`...1`) %>%
#   rename(
#     LD1 = `0`,
#     LD2 = `1`,
#     LD3 = `2`,
#     ground = class
#   )

xform1_confusion <- as.data.frame.matrix(table(xform1$ground, xform1$predicted))
mltools::mcc(confusionM = xform1_confusion)

hand_meas <- read_rds("../data/hand_measurements.rds") %>%
  mutate(
    sub = paste0("sub-", Subject)
  ) %>%
  select(sub, ends_with("_dex"))

# TODO: Rewrite EVR code

# evr <- unlist(ldas$evr) %>%
#   matrix(data = ., nrow = 5, ncol = 3, byrow = TRUE) %>%
#   t() %>%
#   as_tibble() %>%
#   mutate(
#     LD = 1:3
#   )
#
# evr_summary <- evr %>%
#   pivot_longer(-LD, names_to = "group") %>%
#   group_by(LD) %>%
#   summarize(
#     m = 100 * mean(value),
#     sd = 100 * sd(value)
#   ) %>%
#   mutate(
#     se = sd / sqrt(5)
#   )

# Apply LDAs ====

xform <- xform1 %>%
  group_by(hemi) %>%
  mutate(

    LD1z = scale(LD1)[, 1] %>%
      if_else(hemi == "LH", -1 * ., .),

    LD2z = scale(LD2)[, 1],

    LD3z = scale(LD3)[, 1] %>%
      if_else(hemi == "LH", -1 * ., .),

    hand_correct = if_else(str_remove(ground, "..$") == str_remove(predicted, "..$"),
                           "hand_correct", "hand_wrong"),

    label = case_when(
      handedness == "righty" & hand_correct == "hand_correct" ~
        "Righty (correct)",
      handedness == "lefty" & hand_correct == "hand_wrong" ~
        "Lefty (incorrect)",
      handedness == "righty" & hand_correct == "hand_wrong" ~
        "Righty (incorrect)",
      handedness == "lefty" & hand_correct == "hand_correct" ~
        "Lefty (correct)",
    ),

    plot_order = case_when(
      handedness == "righty" & hand_correct == "hand_correct" ~ 1,
      handedness == "lefty" & hand_correct == "hand_wrong" ~ 2,
      handedness == "righty" & hand_correct == "hand_wrong" ~ 3,
      handedness == "lefty" & hand_correct == "hand_correct" ~ 4,
    )

  ) %>%
  arrange(plot_order) %>%
  left_join(hand_meas)

xform_cf <- table(xform$ground, xform$predicted)
xform_error <- sum(diag(xform_cf)) / sum(xform_cf)

xform_chi <- xform %>%
  mutate(
    ground_hand = str_remove(ground, ".H$"),
    ground_hemi = str_extract(ground, ".H"),
    pred_hand = str_remove(predicted, ".H$"),
    pred_hemi = str_extract(predicted, ".H$"),
  ) %>%
  select(ground_hand, ground_hemi, pred_hand, pred_hemi) %>%
  group_by(ground_hand, pred_hand, pred_hemi) %>%
  summarize(
    n = n()
  ) %>%
  pivot_wider(names_from = pred_hand, values_from = n) %>%
  mutate(
    ratio = lefty / (lefty + righty)
  )

lefty_ratio <- (6 + 10) / (6 + 10 + 86 + 82)
righty_ratio <- (17 + 22) / (17 + 22 + 826 + 821)

xform_chi2 <- tibble(true = xform_chi$ratio) %>%
  add_column(expected = c(lefty_ratio, lefty_ratio, righty_ratio, righty_ratio))

chisq.test(xform_chi2, simulate.p.value = TRUE)

xform_sinistral <- xform %>%
  select(sub, handedness, ground, predicted) %>%
  mutate(
    ground_hand = str_remove(ground, ".H$"),
    pred_hand = str_remove(predicted, ".H$"),
    flagged = (ground_hand == "lefty" & pred_hand == "lefty") |
                (ground_hand == "righty" & pred_hand == "lefty")
  ) %>%
  group_by(sub) %>%
  mutate(
    either_flag = sum(flagged) > 0
  ) %>%
  filter(
    sub %in% .$sub[.$either_flag]
  ) %>%
  pivot_wider(id_cols = c(sub, handedness), values_from = pred_hand,
              names_from = hemi)

xform_sinistral %>%
  select(sub, handedness) %>%
  distinct() %>%
  group_by(handedness) %>%
  summarize(
    n = n()
  )

# Summarize for the average icons

xform_summary1 <- xform %>%
  group_by(hemi) %>%
  summarize(
    across(c(LD1, LD2, LD3), list(mean = mean, sd = sd),
           .names = "{.fn}_{.col}")
  )

xform_summary2 <- xform %>%
  group_by(handedness, hemi) %>%
  summarize(
    n = n(),
    across(c(LD1, LD2, LD3), list(mean = mean, sd = sd),
           .names = "{.fn}_{.col}")
  ) %>%
  mutate(
    across(starts_with("sd_"), ~ . / sqrt(n), .names = 'se_{.col}')
  ) %>%
  select_all(~str_replace(., "se_sd_", "se_"))

lefty_color <- "pink2"

## Main publication plot ====


ld1_ld2_plot <- ggplot(mapping = aes(x = LD1, y = LD2)) +
  geom_point(data = xform,
             aes(fill = label, color = label, alpha = label, shape = label),
             size = 2, stroke = 1,
             show.legend = TRUE) +
  scale_shape_manual(values = c(21, 21, 18, 23)) +
  scale_alpha_manual(values = c(1, 1, 0.25, 1)) +
  scale_color_manual(values = c("black","pink3", "grey50", "black")) +
  scale_fill_manual(values = c(lefty_color, lefty_color, "grey50", "grey50")) +
  scale_x_continuous(breaks = c(-5, 0, 5),
                     labels = c("-5 (LH)", "0", "5 (RH)"),
                     sec.axis = dup_axis(breaks = c(-5, 5),
                                         labels = c(expression(Z(LD1)[LH]*" = 0"),
                                                    expression(Z(LD1)[RH]*" = 0")),
                                         name = NULL)) +
  theme_bw() +
  labs(x = "LD1 (95%)", y = "LD2 (3%)",
       color = "Handedness-1", fill = "Handedness-1", shape = "Handedness-1",
       alpha = "Handedness-1",
       # caption = "All hemispheres were classified correctly."
  ) +
  theme(legend.position = "bottom") +
  guides(
    shape = guide_legend(override.aes = list(alpha = 1),
                         nrow = 2, byrow = FALSE)
  ) +
  annotate(geom = "point",
           x = xform_summary2$mean_LD1, y = xform_summary2$mean_LD2,
           shape = 22, size = 4.5,
           fill = c(lefty_color, lefty_color, "grey50", "grey50"),
           stroke = 1
  ) +
  annotate(geom = "segment",
           x = xform_summary2$mean_LD1,
           y = xform_summary2$mean_LD2 - xform_summary2$se_LD2,
           yend = xform_summary2$mean_LD2 + xform_summary2$se_LD2) +
  annotate(geom = "segment",
           y = xform_summary2$mean_LD2,
           x = xform_summary2$mean_LD1 - xform_summary2$se_LD1,
           xend = xform_summary2$mean_LD1 + xform_summary2$se_LD1)

ggsave("plots/ld1_ld2_size-manuscript.png", width = 6, height = 4, units = "in")

ld1_ld2_plot_poster <- ggplot(xform, mapping = aes(x = LD1, y = LD2)) +
  geom_point(aes(color = handedness, alpha = handedness),
             size = 3, stroke = 1,
             show.legend = TRUE) +
  scale_x_continuous(breaks = c(-5, 0, 5),
                     labels = c("-5 (LH)", "0", "5 (RH)")) +
  scale_color_manual(values = c("violet", "black")) +
  scale_alpha_manual(values = c(0.8, 0.2)) +
  theme_bw() +
  labs(x = "LD1 (95%)", y = "LD2 (3%)",
       color = "Handedness-1", fill = "Handedness-1", shape = "Handedness-1",
       alpha = "Handedness-1",
       # caption = "All hemispheres were classified correctly."
  ) +
  annotate(geom = "point",
           x = xform_summary2$mean_LD1, y = xform_summary2$mean_LD2,
           shape = 22, size = 6,
           fill = c('violet', "violet", "grey50", "grey50"),
           stroke = 2
  ) +
  theme(legend.position = "bottom",
        text = element_text(size = 24))


ggsave(ld1_ld2_plot_poster, filename = "plots/ld1_ld2_size-poster.png",
       width = 10, height = 8, units = "in")


# Simple plot for research statement ====


## LD2 and LD3 ====

lm(LD2 ~ EHI, data = xform) %>%
  summary()

lm(LD2z ~ EHI, data = xform) %>%
  summary()

lm_LD3z <- lm(LD3z ~ EHI, data = xform)

ggplot(xform, aes(x = EHI, y = LD3z, color = hemi)) +
  scale_size() +
  geom_smooth(method = "lm") +
  theme_bw()


## Z(LD1) ====

# Descriptive labels for interpretation of re-scaled values
desc_labels <- c("Closer to\ncontralateral", "Average", "Super\nipsilateral")

ggplot(xform, aes(x = jitter(EHI), y = LD1z, color = hemi)) +
  geom_point(alpha = 0.25, shape = 18) +
  geom_vline(xintercept = c(0, 30), color = "red") +
  geom_smooth(method = "lm") +
  scale_y_continuous(breaks = seq(-4, 4, by = 2),
                     sec.axis = sec_axis(~ .,
                                         breaks = c(-2.5, 0, 2.5),
                                         labels = desc_labels)) +

  theme_bw() +
  labs(x = "EHI (jittered)", y = "Z(LD1)",
       title = "Linear regression")

ggplot(xform, aes(x = jitter(EHI), y = LD1z, color = hemi)) +
  geom_point(alpha = 0.25, shape = 18) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2)) +
  scale_y_continuous(breaks = seq(-4, 4, by = 2),
                     sec.axis = sec_axis(~ .,
                                         breaks = c(-2.5, 0, 2.5),
                                         labels = desc_labels)) +

  theme_bw() +
  labs(x = "EHI (jittered)", y = "Z(LD1)",
       title = "Polynomial regression")

ggplot(xform, aes(x = LI_dex, y = LD1z, color = hemi)) +
  geom_point(alpha = 0.25, shape = 18) +
  geom_smooth(method = "lm") +
  scale_y_continuous(breaks = seq(-4, 4, by = 2),
                     sec.axis = sec_axis(~ .,
                                         breaks = c(-2.5, 0, 2.5),
                                         labels = desc_labels)) +
  theme_bw() +
  labs(x = "Dexterity LI", y = "Z(LD1)")

lm(LD1z ~ LI_dex*hemi, data = xform) %>%
  summary()

xform_long <- xform %>%
  pivot_longer(matches("^LD.$"), names_to = "LD")

ggplot(xform_long, aes(y = value, fill = gender)) +
  geom_boxplot() +
  facet_wrap(vars(LD), scales = "free") +
  theme_bw()

## Z(LD1) for both hemispheres ====

xform2 <- xform %>%
  select(sub, gender, handedness, handedness2, EHI, hemi, LD1z) %>%
  pivot_wider(names_from = hemi, values_from = LD1z) %>%
  mutate(
    plot_order = if_else(handedness == "righty", 1, 2)
  ) %>%
  arrange(plot_order)

xform2_cor <- xform2 %>%
  group_by(handedness) %>%
  nest() %>%
  mutate(
    r = map(data, ~cor.test(.x$LH, .x$RH))
  )

xform2_cor2 <- xform2 %>%
  group_by(handedness2) %>%
  nest() %>%
  mutate(
    r = map(data, ~cor.test(.x$LH, .x$RH))
  )


png("plots/LD1_LHvRH.png", width = 5, height = 4, res = 300, units = "in")

ggplot(xform2, aes(x = LH, y = RH)) +
  geom_point(aes(color = handedness), alpha = 0.5) +
  geom_smooth(aes(color = handedness), method = "lm") +
  theme_bw() +
  labs(x = "LH LD1 value", y = "RH LD1 value",
       title = "Correlation between LD1 values within-participant")

dev.off()

lh_rh_lm <- xform2 %>%
  select(LH, RH, handedness) %>%
  mutate(
    lefty = handedness == "lefty"
  ) %>%
  lm(RH ~ LH*lefty, data = .)

summary(lh_rh_lm)
coef(lh_rh_lm)[2] + coef(lh_rh_lm)[4]

ggplot(xform2, aes(x = LH, y = RH)) +
  geom_point(aes(color = gender), alpha = 0.5) +
  geom_smooth(aes(color = gender), method = "lm") +
  theme_bw()

## Segmented linear regression ====

library(segmented)

# Does EHI predict Z(LD1) score?
lm_intonly <- lm(LD1z ~ 1, data = xform)
lm_ehi <- lm(LD1z ~ EHI, data = xform)

anova(lm_intonly, lm_ehi) # YES

# Does hemisphere affect LD1z score?
lm_ehihemi <- lm(LD1z ~ EHI*hemi, data = xform)
anova(lm_ehi, lm_ehihemi) # NO

# Does gender affect this?
lm_gender <- lm(LD1z ~ EHI*gender, data = xform)
anova(lm_ehi, lm_gender) # NO

## Segmented ====

lm_segn1 <- segmented::segmented(lm_ehi, ~EHI, npsi = 1)
lm_segn2 <- segmented::segmented(lm_ehi, ~EHI, npsi = 2)

anova(lm_ehi, lm_segn1, lm_segn2)

AIC(lm_ehi, lm_segn1, lm_segn2)

rmses <- sapply(list(lm_intonly, lm0, lm1, lm1_poly, lm2),
                function(x) sqrt(mean(x$residuals^2)))

lm_LD3z_segn1 <- segmented::segmented(lm_LD3z, ~EHI, npsi = 1)
lm_LD3z_segn2 <- segmented::segmented(lm_LD3z, ~EHI, npsi = 2)
anova(lm_LD3z, lm_LD3z_segn1, lm_LD3z_segn2)

## Quadratic ====

lm_quad <- lm(LD1z ~ EHI + I(EHI^2), data = xform)

anova(lm_ehi, lm_quad)



# Plots

png("plots/ZLD1_segn1.png", width = 6, height = 4, res = 300, units = "in")

p1 <- ggplot(xform, aes(x = jitter(EHI), y = LD1z)) +
  geom_point(aes(color = hemi), alpha = 0.5, shape = 18, size = 2) +
  scale_color_manual(values = c("cyan4", "brown2")) +
  geom_vline(xintercept = c(75, 75 + c(-7.432, 7.432)), color = "black",
             linetype = c("solid", "dotted", "dotted")) +
  geom_smooth(
    data = filter(xform, EHI <= 75),
    method = "lm"
  ) +
  geom_smooth(
    data = filter(xform, EHI > 75),
    method = "lm"
  ) +
  theme_bw() +
  labs(x = "EHI (jittered)", y = "Z(LD1)", color = "Hemisphere") +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 4)))

ggMarginal(p1, type = "histogram")

dev.off()

png("plots/ZLD3_segn1.png", width = 6, height = 4, res = 300, units = "in")

p_LD3z <- ggplot(xform, aes(x = jitter(EHI), y = LD3z)) +
  geom_point(aes(color = hemi), alpha = 0.5, shape = 18, size = 2) +
  scale_color_manual(values = c("cyan4", "brown2")) +
  geom_vline(xintercept = c(75, 75 + c(-8.042, 8.042)), color = "black",
             linetype = c("solid", "dotted", "dotted")) +
  geom_smooth(
    data = filter(xform, EHI <= 75),
    method = "lm"
  ) +
  geom_smooth(
    data = filter(xform, EHI > 75),
    method = "lm"
  ) +
  theme_bw() +
  labs(x = "EHI (jittered)", y = "Z(LD3)", color = "Hemisphere") +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 4)))

ggMarginal(p_LD3z, type = "histogram")

dev.off()

## Research statement plots

ggplot(xform, aes(x = jitter(EHI), y = LD1z)) +
  geom_point(aes(color = hemi), alpha = 0.5, shape = 18, size = 2) +
  geom_smooth(
    data = filter(xform, EHI <= 75),
    method = "lm"
  ) +
  geom_smooth(
    data = filter(xform, EHI > 75),
    method = "lm"
  ) +
  theme_bw() +
  labs(x = "EHI (jittered)", y = "Z(LD1)", color = "Hemisphere") +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 4)))

## Quickly pull EHI values in order to get descriptive statistics ====
EHI <- hc %>%
  select(sub, gender, EHI) %>%
  distinct() %>%
  left_join(hand_meas, join_by("sub")) %>%
  mutate(
    dominant_dex = if_else(EHI > 0, Right_dex, Left_dex),
    nondom_dex = if_else(EHI > 0, Left_dex, Right_dex)
  )

EHI %>%
  group_by(gender) %>%
  summarize(
    n = n(),
    across(c(EHI, dominant_dex, nondom_dex),
           list(\(x) mean(x, na.rm = TRUE), \(x) sd(x, na.rm = TRUE))
           )
  )

median(EHI$EHI)
mean(EHI$EHI)

table(EHI$EHI > 0)
table(EHI$EHI > 30)

lm(EHI ~ gender, data = EHI) %>% summary()

lm(Right_dex ~ gender, data = EHI) %>% summary()
lm(Left_dex ~ gender, data = EHI) %>% summary()
lm(dominant_dex ~ gender, data = EHI) %>% summary()
lm(nondom_dex ~ gender, data = EHI) %>% summary()

lm(LI_dex ~ gender, data = EHI) %>% summary()



ggplot(data = EHI, aes(x = EHI, y = LI_dex, color = gender)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw()

# Hemisphere difference ====

wide <- xform %>%
  select(sub, gender, LD1, LD2, LD3, EHI, handedness, hemi) %>%
  pivot_wider(names_from = hemi, values_from = c(LD1, LD2, LD3)) %>%
  mutate(
    euclidean = sqrt((LD1_RH - LD1_LH)^2 + (LD2_RH - LD1_LH)^2 +
                       (LD3_RH - LD3_LH)^2),
    LD1_dist = abs(LD1_LH) + LD1_RH
  )

wide_lm <- lm(euclidean ~ LD1_dist, data = wide)
wide_lm_gender <- lm(euclidean ~ LD1_dist + gender, data = wide)

anova(wide_lm, wide_lm_gender)

wide_aug <- broom::augment(wide_lm)

# Are lefties's eucliean distance farther apart than LD1 would suggest?
wide <- left_join(wide, select(wide_aug, LD1_dist, euclidean, .resid))

cor(wide$euclidean, wide$LD1_dist) # r=.95

ggplot(wide, aes(x = EHI, y = euclidean)) +
  geom_point() +
  geom_smooth() +
  theme_bw() +
  labs(title = "Euclidean (LD1,2,3) distance)")

ggplot(wide, aes(x = EHI, y = LD1_dist)) +
  geom_point() +
  geom_smooth() +
  theme_bw() +
  labs(title = "Horizontal distance")

## Research statement ====

random_sub <- xform %>%
  ungroup() %>%
  filter(
    sub %in% sample(xform$sub, 1)
    # sub == "sub-105620"
  )

half_dist <- abs(diff(random_sub$LD1)) / 2
center <- min(random_sub$LD1) + half_dist

res_statement1 <- (
    ggplot(xform, aes(x = LD1, y = sub)) +
      geom_point(aes(color = hemi),
                 size = 1.5, alpha = 0.5, shape = 20,
                 show.legend = TRUE) +
      geom_line(
        data = random_sub, linetype = "solid"
      ) +
      geom_label(x = center, y = unique(random_sub$sub), label = "dist.") +
      scale_x_continuous(breaks = c(-5, 0, 5),
                         labels = c("-5 (LH)", "0", "5 (RH)")) +
      scale_y_discrete(labels = NULL) +
      theme_classic() +
      theme(text = element_text(size = 10.5),
            legend.position = "none",
            axis.ticks.y = element_blank()) +
      labs(x = "Hemisphere similarity", y = "Participant", color = "Hemi.")
  )


res_statement2 <- ggplot(wide, aes(x = jitter(EHI), y = LD1_dist)) +
  geom_point(alpha = 0.25, shape = 20) +
  geom_smooth(se = FALSE) +
  scale_x_continuous(breaks = c(-100, 0, 100),
                     labels = c("-100 (left)", "0 (ambi.)", "100 (right)")) +
  scale_y_continuous(breaks = c(4, 8, 12),
                     labels = c("Similar", "", "Dissimilar")) +
  theme_classic() +
  theme(text = element_text(size = 10.5)) +
  labs(x = "Handedness score", y = "Distance")

library(patchwork)

png("res_statement.png", width = 5.5, height = 1.75, units = "in",
    res = 300)

res_statement1 + res_statement2

dev.off()

## Models ====

ld1_dist_lm_intonly <- lm(LD1_dist ~ 1, data = wide)
ld1_dist_lm_ehi <- lm(LD1_dist ~ EHI, data = wide)

# No effect of gender
anova(ld1_dist_lm_intonly, ld1_dist_lm_ehi)

ld1_dist_lm_segn1 <- segmented::segmented(ld1_dist_lm_ehi, ~EHI, npsi = 1)
ld1_dist_lm_segn2 <- segmented::segmented(ld1_dist_lm_ehi, ~EHI, npsi = 2)

anova(ld1_dist_lm_intonly, ld1_dist_lm_segn1, ld1_dist_lm_segn2)

AIC(ld1_dist_lm_intonly, ld1_dist_lm_segn1, ld1_dist_lm_segn2)

ld1_dist_lm_quad <- lm(LD1_dist ~ EHI + I(EHI^2), data = wide)

anova(ld1_dist_lm_intonly, ld1_dist_lm_quad)
AIC(ld1_dist_lm_intonly, ld1_dist_lm_quad)

quad_peak_x <- -coef(ld1_dist_lm_quad)[2] / (2 * coef(ld1_dist_lm_quad)[3])
quad_peak_y <- predict(ld1_dist_lm_quad, tibble(EHI = quad_peak))


ggplot(wide, aes(x = jitter(EHI), y = LD1_dist)) +
  geom_point(alpha = 0.5) +
  geom_smooth(
    method = "lm",
    formula = y ~ poly(x, 2)
  ) +
  annotate("point", x = quad_peak, y = quad_peak_y, color = "red", size = 2.5) +
  theme_bw() +
  labs(x = "EHI (jittered)", y = "|LH| - RH",
       title = "Horizontal distance")

ggsave("plots/wide_quadratic.png", width = 6, height = 4, units = "in")

ggplot(wide, aes(x = euclidean, y = LD1_dist)) +
  geom_point(aes(fill = EHI), shape = 21) +
  viridis::scale_fill_viridis() +
  geom_smooth(aes(color = handedness)) +
  theme_bw()

ggplot(wide, aes(x = EHI, y = .resid)) +
  geom_point(aes(fill = EHI), shape = 21) +
  viridis::scale_fill_viridis() +
  geom_smooth(method = "lm") +
  theme_bw()

png("plots/difference_LD1.png", width = 6, height = 4, res = 300, units = "in")

ggplot(wide, aes(x = jitter(EHI), y = LD1_dist)) +
  geom_point(alpha = 0.5) +
  geom_smooth(
    data = filter(wide, EHI <= 75),
    method = "lm"
  ) +
  geom_smooth(
    data = filter(wide, EHI > 75),
    method = "lm"
  ) +
  theme_bw() +
  labs(x = "EHI (jittered)", y = "|LH| - RH",
       title = "Horizontal distance")

dev.off()

# Aux plots ====

## LD1, LD3 ====

png("plots/ld1_ld3.png", width = 6, height = 4, res = 300, units = "in")

ggplot(mapping = aes(x = LD1, y = LD3)) +
  geom_point(data = xform,
             aes(fill = label, color = label, alpha = label, shape = label),
             size = 2, stroke = 1,
             show.legend = TRUE) +
  scale_shape_manual(values = c(21, 21, 18, 23)) +
  scale_alpha_manual(values = c(1, 1, 0.25, 1)) +
  scale_color_manual(values = c("black","pink3", "grey50", "black")) +
  scale_fill_manual(values = c(lefty_color, lefty_color, "grey50", "grey50")) +
  scale_x_continuous(breaks = c(-5, 0, 5),
                     labels = c("-5 (LH)", "0", "5 (RH)")) +
  theme_bw() +
  labs(x = "LD1 (95%)", y = "LD3 (2%)",
       color = "Handedness-1", fill = "Handedness-1", shape = "Handedness-1",
       alpha = "Handedness-1",
       # caption = "All hemispheres were classified correctly."
  ) +
  theme(legend.position = "bottom") +
  guides(
    shape = guide_legend(override.aes = list(alpha = 1),
                         nrow = 2, byrow = FALSE)
  ) +
  annotate(geom = "point",
           x = xform_summary2$mean_LD1, y = xform_summary2$mean_LD3,
           shape = 22, size = 4.5,
           fill = c(lefty_color, lefty_color, "grey50", "grey50"),
           stroke = 1
  ) +
  annotate(geom = "segment",
           x = xform_summary2$mean_LD1,
           y = xform_summary2$mean_LD3 - xform_summary2$se_LD3,
           yend = xform_summary2$mean_LD3 + xform_summary2$se_LD3) +
  annotate(geom = "segment",
           y = xform_summary2$mean_LD3,
           x = xform_summary2$mean_LD1 - xform_summary2$se_LD1,
           xend = xform_summary2$mean_LD1 + xform_summary2$se_LD1)

dev.off()

# LD3 plots ====

png("plots/ld2_ld3.png", width = 4, height = 5, units = "in", res = 300)

ggplot(mapping = aes(x = LD2, y = LD3)) +
  geom_point(data = xform,
             aes(fill = label, color = label, alpha = label, shape = label),
             size = 2, stroke = 1,
             show.legend = TRUE) +
  scale_shape_manual(values = c(21, 21, 18, 23)) +
  scale_alpha_manual(values = c(1, 1, 0.25, 1)) +
  scale_color_manual(values = c("black","pink3", "grey50", "black")) +
  scale_fill_manual(values = c(lefty_color, lefty_color, "grey50", "grey50")) +
  theme_bw() +
  labs(x = "LD2 (3%)", y = "LD3 (2%)",
       color = "Handedness-1", fill = "Handedness-1", shape = "Handedness-1",
       alpha = "Handedness-1",
       # caption = "All hemispheres were classified correctly."
  ) +
  theme(legend.position = "bottom") +
  guides(
    shape = guide_legend(override.aes = list(alpha = 1),
                         nrow = 2, byrow = FALSE)
  )

dev.off()
