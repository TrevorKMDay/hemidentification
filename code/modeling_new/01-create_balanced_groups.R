library(tidyverse)
library(groupdata2)

library(MatchIt)

library(patchwork)

select <- dplyr::select

setwd(str_glue("~/Google Drive/My Drive/Projects/hemisphere_fingerprinting/",
               "code/modeling_new/"))

# Demographic data ====

sex <- read_csv("../../data/unrestricted_data.csv", show_col_types = FALSE) %>%
  mutate(
    sub = paste0("sub-", Subject),
  ) %>%
  select(sub, Gender, Age) %>%
  rename(
    gender = Gender,
    age_group = Age
  )

hands <- read_csv("../../data/restricted_data.csv", show_col_types = FALSE) %>%
  rename(
    EHI = Handedness,
    age = Age_in_Yrs,
    family = Family_ID
  ) %>%
  mutate(
    sub = paste0("sub-", Subject),
    handedness = if_else(EHI > 0, "righty", "lefty"),
  ) %>%
  select(sub, family, EHI, age, starts_with("handedness"))

# Merge data
demo <- left_join(sex, hands, by = join_by(sub))

ggplot(demo, aes(x = EHI, fill = gender)) +
  geom_histogram(binwidth = 10, boundary = 5) +
  scale_x_continuous(breaks = seq(-100, 100, 20)) +
  theme_bw()

ggplot(demo, aes(x = EHI, fill = gender)) +
  geom_boxplot() +
  scale_x_continuous(breaks = seq(-100, 100, 20)) +
  theme_bw()


# Get the list of subjects to keep =====

pconns <- list.files("../../data/pconns_0311_complete/", full.names = TRUE)

pconns_language <- list.files("../../data/pconns_0424_language/",
                              full.names = TRUE)

subs_to_keep_rest <- unique(str_remove(basename(pconns), "_.*.txt"))
subs_to_keep_lang <- unique(str_remove(basename(pconns_language), "_.*.txt"))

subs_in_either <- unique(c(subs_to_keep_rest, subs_to_keep_lang))

sum(subs_to_keep_rest %in% subs_to_keep_lang)

# Create groups based on resting-state data ====

## These groups are carried forward into all manipulations

# This does two folds, one group fold balanced on gender/age group,
# because it can only handle one categorical grouping variable and
# continuous EHI because EHI is more important to our question.

# The second fold is which hemisphere to keep if only one hemisphere
# is being kept.

set.seed(20007)

input <- demo %>%
  filter(
    sub %in% subs_in_either,
    # The very small 36+ group when divided into
    # age groups didn't work well
    # age_group != "36+"
  ) %>%
  mutate(
    sub = factor(sub),
    family = factor(family),
    # gender_age = factor(paste(gender, age_group))
  ) 

# Set up folding 

task_folds1 <- as_task_classif(input, target = "EHI")
task_folds1$col_roles$group <- "family"
task_folds1$col_roles$stratum <- c("gender", "age_group")

cv5 <- rsmp("same_other_sizes_cv", folds = 5)
set.seed(20007)
cv5$instantiate(task_folds1)
five_fold <- cv5$instance

cv2 <- rsmp("same_other_sizes_cv", folds = 2)
set.seed(20007)
cv2$instantiate(task_folds1)
hemi_fold <- cv2$instance

groups <- input %>%
  add_column(
    fold = five_fold$fold.dt$fold,
    hemi_fold = hemi_fold$fold.dt$fold
  ) %>%
  mutate(

    # Convert to A-E
    group = as.character(fold) %>%
      recode_values(from = as.character(1:5), to = LETTERS[1:5]) %>%
      factor(),

    # Convert to LH/RH
    hemi_fold = factor(recode_values(hemi_fold, from = 1:2, to = c("LH", "RH")))

  ) %>%
  select(-fold)

# Check families

families <- groups %>%
  group_by(family) %>%
  nest() %>%
  mutate(
    n = map_int(data, nrow),
    groups = map_int(data, ~length(unique(.x$group)))
  )

groups_summary <- groups %>%
  group_by(group) %>%
  summarize(
    n = n(),
    across(c(age, EHI), c("mean" = mean, "sd" = sd)),
    propM_mean = sum(gender == "M") / n(),
  ) %>%
  mutate(
    propM_se = sqrt((propM_mean) * (propM_mean) / n),
    across(ends_with("sd"), ~ .x / sqrt(n), .names = "{.col}_se")
  )

## Break down by gender ====

bygender_summary <- groups %>%
  group_by(group, gender) %>%
  summarize(
    n = n(),
    across(c(age, EHI), c("mean" = mean, "sd" = sd)),
  ) %>%
  mutate(
    across(ends_with("sd"), ~ .x / sqrt(n), .names = "{.col}_se")
  )

p1 <- ggplot(bygender_summary, aes(x = group, y = age_mean, color = gender)) +
  geom_pointrange(aes(ymin = age_mean - age_sd_se,
                      ymax = age_mean + age_sd_se),
                  position = position_dodge2(width = .25),
                  size = 0.25) +
  scale_y_continuous(limits = c(22, 35), breaks = seq(22, 35, by = 2)) +
  theme_bw() +
  labs(x = "Fold", y = "Age in years (SE)", title = "(A)", color = "Gender") +
  theme(legend.position = "bottom")

p2 <- ggplot(bygender_summary, aes(x = group, y = EHI_mean, color = gender)) +
  geom_pointrange(aes(ymin = EHI_mean - EHI_sd_se,
                      ymax = EHI_mean + EHI_sd_se),
                  position = position_dodge2(width = .25),
                  size = 0.25) +
  scale_y_continuous(limits = c(-100, 100)) +
  theme_bw() +
  labs(x = "Fold", y = "EHI (SE)", title = "(B)", color = "Gender") +
  theme(legend.position = "bottom")

## Breakdown by handedness group ====

byhand_summary <- groups %>%
  mutate(
    righthanded = case_when(
      EHI > 0 ~ "right",
      EHI <= 0 ~ "left",
    )
  ) %>%
  group_by(group, righthanded) %>%
  summarize(
    n = n(),
    propM_mean = 100 * sum(gender == "M") / n(),
  ) %>%
  mutate(
    propM_se = sqrt((propM_mean) * (propM_mean) / n),
  )

p3 <- ggplot(byhand_summary, aes(x = group, y = propM_mean,
                                 color = righthanded)) +
  geom_pointrange(aes(ymin = propM_mean - propM_se,
                      ymax = propM_mean + propM_se,),
                  position = position_dodge2(width = .25),
                  size = 0.25) +
  scale_y_continuous(limits = c(0, 100)) +
  scale_color_manual(values = c("orange", "purple")) +
  theme_bw() +
  labs(x = "Fold", y = "% male (SE)", title = "(C)",
       color = "Handedness") +
  theme(legend.position = "bottom")


p1 + p2 + p3 +
  plot_layout(guides = "collect") &
  theme(legend.position = 'bottom')

ggsave("plots/demographics1.png", width = 6.5, height = 4)

t.test(age ~ gender, data = demo)

## Hemi check ====

bygender_byhemi_summary <- groups %>%
  group_by(hemi_fold, gender) %>%
  summarize(
    n = n(),
    across(c(age, EHI), c("mean" = mean, "sd" = sd)),
  ) %>%
  mutate(
    across(ends_with("sd"), ~ .x / sqrt(n), .names = "{.col}_se")
  )

p1h <- ggplot(bygender_byhemi_summary,
             aes(x = hemi_fold, y = age_mean, color = gender)) +
  geom_pointrange(aes(ymin = age_mean - age_sd_se,
                      ymax = age_mean + age_sd_se),
                  position = position_dodge2(width = .25),
                  size = 0.25) +
  scale_y_continuous(limits = c(22, 35), breaks = seq(22, 35, by = 2)) +
  theme_bw() +
  labs(x = "Fold", y = "Age in years (SE)", title = "(A)", color = "Gender") +
  theme(legend.position = "bottom")

p2h <- ggplot(bygender_byhemi_summary,
              aes(x = hemi_fold, y = EHI_mean, color = gender)) +
  geom_pointrange(aes(ymin = EHI_mean - EHI_sd_se,
                      ymax = EHI_mean + EHI_sd_se),
                  position = position_dodge2(width = .25),
                  size = 0.25) +
  scale_y_continuous(limits = c(-100, 100)) +
  theme_bw() +
  labs(x = "Fold", y = "EHI (SE)", title = "(B)", color = "Gender") +
  theme(legend.position = "bottom")

byhand_byhemi_summary <- groups %>%
  mutate(
    righthanded = case_when(
      EHI > 0 ~ "right",
      EHI <= 0 ~ "left",
    )
  ) %>%
  group_by(hemi_fold, righthanded) %>%
  summarize(
    n = n(),
    propM_mean = 100 * sum(gender == "M") / n(),
  ) %>%
  mutate(
    propM_se = sqrt((propM_mean) * (propM_mean) / n),
  )

p3h <- ggplot(byhand_byhemi_summary,
             aes(x = hemi_fold, y = propM_mean, color = righthanded)) +
  geom_pointrange(aes(ymin = propM_mean - propM_se,
                      ymax = propM_mean + propM_se,),
                  position = position_dodge2(width = .25),
                  size = 0.25) +
  scale_y_continuous(limits = c(0, 100)) +
  scale_color_manual(values = c("orange", "purple")) +
  theme_bw() +
  labs(x = "Fold", y = "% male (SE)", title = "(C)",
       color = "Handedness") +
  theme(legend.position = "bottom")

p1h + p2h + p3h +
  plot_layout(guides = "collect") &
  theme(legend.position = 'bottom')

ggsave("plots/demographics2.png", width = 6.5, height = 4)

# Write out  ====

write_csv(select(groups, -family), "folds_to_use2.csv")

# Run models to show balance ====

library(multcomp)

## Folds =====

age_anova <- aov(age ~ group, data = groups)
summary(age_anova)
# apa_print(age_anova)

ehi_test <- kruskal.test(EHI ~ group, data = groups)

gender_test <- glm(gender == "M" ~ group, family = binomial, data = groups) 
summary(glht(gender_test, mcp(group = "Tukey")))

## Hemi folds ====

age_anova_hemi <- aov(age ~ hemi_fold, data = groups)
summary(age_anova_hemi)
# apa_print(age_anova_hemi)

ehi_test_hemi <- kruskal.test(EHI ~ hemi_fold, data = groups)

gender_test_hemi <- glm(gender == "M" ~ hemi_fold, family = binomial, 
                        data = groups) 

summary(gender_test_hemi)

# Groups for later =====

right_handed <- groups %>%
  filter(
    EHI > 0
  )

sum(right_handed$gender == "M") / nrow(right_handed)
mean(right_handed$age)

left_handed <- groups %>%
  filter(
    EHI <= 0
  )

sum(left_handed$gender == "M") / nrow(left_handed)
mean(left_handed$age)

left_handed_groups <- left_handed %>%
  group_by(group) %>%
  nest() %>%
  mutate(
    n = map_int(data, nrow)
  )

set.seed(20007)

# Randomly match groups 
right_handed_random <- right_handed %>%
  group_by(group) %>%
  nest() %>%
  left_join(
    select(left_handed_groups, group, n)
  ) %>%
  mutate(
    subsample = map(data, ~sample_n(.x, size = n))
  ) %>%
  select(-data) %>%
  unnest(subsample)




left_handed2 <- left_handed %>%
  rowwise() %>%
  mutate(
    dummytreat = sample(0:1, 1)
  )

match_me <- groups %>%
  mutate(
    absEHI = abs(EHI),
    lefty = EHI <= 0
  )

m.out <- matchit(lefty ~ gender + age + absEHI, 
                 data = match_me,
                 method = "nearest")

m.data <- match.data(m.out) %>%
  arrange(subclass)

right_handed_matched1 <- groups %>%
  mutate(
    absEHI = abs(EHI),
    lefty = EHI <= 0
  ) %>%
  group_by(group) %>%
  nest() %>%
  mutate(
    match_obj = map(data, ~matchit(lefty ~ gender + age + absEHI, 
                                    data = .x,
                                    method = "nearest")),
    match_data = map(match_obj, match.data)
  )

right_handed_matched2 <- right_handed_matched1 %>%
  select(group, match_data) %>%
  unnest(match_data) %>%
  filter(
    # !lefty 
  ) 

plot_me <- right_handed_matched2 %>%
  mutate(
    lefty = if_else(lefty, "lefty", "righty_matched")
  ) %>%
  bind_rows(
    right_handed_random %>%
      filter(EHI > 0) %>%
      mutate(lefty = "random")
  )

plot_me %>%
  group_by(lefty) %>%
  summarize(
    n = n(),
    var = var(abs(EHI))
  )

ggplot(plot_me, aes(x = lefty, y = abs(EHI), fill = lefty)) +
  geom_boxplot() +
  theme_bw()

# Write out simple table
 
righty_matched_scored <- right_handed_matched2 %>%
  filter(
    EHI > 0
  ) %>%
  select(sub, group) %>%
  mutate(
    paired = TRUE
  )

righty_matched_groups <- right_handed_random %>%
  select(sub, group) %>%
  mutate(
    random = TRUE
  ) 

righty_matched_groups_final <- full_join(righty_matched_scored, 
                                          righty_matched_groups) %>%
  arrange(sub) %>%
  mutate(
    across(everything(), ~replace_na(.x, FALSE))
  )

write_csv(righty_matched_groups_final, "righty_subsample_subs.csv")

# Plot matched groups

matching_summary <- plot_me %>%
  group_by(lefty, group) %>%
  summarize(
    n = n(),
    n_male = sum(gender == "M"),
    mean_age = mean(age),
    sd_age = sd(age),
    mean_EHI = mean(abs(EHI)),
    sd_EHI = sd(abs(EHI))
  ) %>%
  mutate(
    pct_male = n_male / n,
    se_pct_male = sqrt((pct_male * (1 - pct_male)) / n),
    se_age = sd_age / sqrt(n),
    se_EHI = sd_EHI / sqrt(n),

    lefty = recode_values(
      lefty, 
      "lefty" ~ "Sinistrals",
      "random" ~ "Dextral-random",
      "righty_matched" ~ "Dextral-matched"
    )

  )

matching_age <- ggplot(matching_summary, aes(x = group, color = lefty)) +
  geom_pointrange(aes(y = mean_age, ymin = mean_age - se_age, 
                      ymax = mean_age + se_age),
                  position = position_dodge2(width = 0.25)) +
  scale_y_continuous(limits = c(22, 34)) +
  theme_bw() +
  labs(x = "Fold", y = "Age in years (SE)", title = "(A)", color = "Group")

matching_ehi <- ggplot(matching_summary, aes(x = group, color = lefty)) +
  geom_pointrange(aes(y = mean_EHI, ymin = mean_EHI - se_EHI, 
                      ymax = mean_EHI + se_EHI),
                  position = position_dodge2(width = 0.25)) +
  scale_y_continuous(limits = c(0, 100)) +
  theme_bw() +
  labs(x = "Fold", y = "Abs. EHI (SE)", title = "(B)", color = "Group")

matching_gender <- ggplot(matching_summary, aes(x = group, color = lefty)) +
  geom_pointrange(aes(y = pct_male, ymin = pct_male - se_pct_male, 
                      ymax = pct_male + se_pct_male),
                  position = position_dodge2(width = 0.25)) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw() +
  labs(x = "Fold", y = "% male (SE)", title = "(C)", color = "Group")

matching_age + matching_ehi + matching_gender + 
  plot_layout(guides = "collect") &
  theme(legend.position = 'bottom')

ggsave("plots/demographics3_subsampledrighties.png", width = 6.5, height = 4)

# Balanced model tests

age_anova_blancing <- aov(age ~ group + lefty, data = plot_me)
summary(age_anova_blancing)
papaja::apa_print(age_anova_blancing)

plot_me2 <- plot_me %>%
  mutate(
    absEHI = abs(EHI),
    group2 = paste0(lefty, group)
  )

# ehi_test_balancing <- kruskal.test(absEHI ~ group2, data = plot_me2)

ehi_test_balancing_SRH <- rcompanion::scheirerRayHare(absEHI ~ lefty + group, 
                                                      data = plot_me2)

gender_test_balancing <- glm(gender == "M" ~ lefty + group, family = binomial, 
                            data = plot_me) 

summary(gender_test_balancing)
