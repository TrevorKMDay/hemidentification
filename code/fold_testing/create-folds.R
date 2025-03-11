library(tidyverse)
library(here)
library(randomizr)

# Load data ====

dat0 <- read_csv(here("data", "unrestricted_data.csv"),
                 show_col_types = FALSE) %>%
  select_all(~str_remove_all(., "-|3T_"))

rdat0 <- read_csv(here("data", "restricted_data.csv"), show_col_types = FALSE)

## Tidy up data ====

dat <- dat0 %>%
  mutate(
    anat_include = if_else(T1_Count > 0 & T2_Count > 0, "anat_include",
                           "anat_exclude"),
    rs_include = if_else(RSfMRI_PctCompl >= 50, "rs_include", "rs_exclude"),
    task_include = if_else(Full_Task_fMRI, "task_include", "task_exclude"),
    rstask_include = if_else(rs_include == "rs_include" &
                              task_include == "task_include",
                            "rstask_include", "rstask_exclude")
  ) %>%
  select(Subject, Gender, matches("^T._Count"), RSfMRI_PctCompl,
         Full_Task_fMRI, ends_with("include"))

ggplot(dat, aes(x = RSfMRI_PctCompl)) +
  geom_histogram(boundary = 0, binwidth = 5) +
  theme_bw()

dat_exclusions <- dat %>%
  select(Subject, Gender, ends_with("include")) %>%
  rename(
    sub = Subject,
    gender = Gender
  )

rdat <- rdat0 %>%
  select(Subject, Age_in_Yrs, Family_ID, Handedness) %>%
  rename(
    sub = Subject,
    age = Age_in_Yrs,
    family = Family_ID,
    ehi = Handedness
  ) %>%
  mutate(
    hand1 = if_else(ehi > 0, "right", "left"),
    hand2 = if_else(ehi > 50, "right2", "non-right")
  )

# Everyone with sufficient anatomical and RS data, butpossibly insufficient
#   task data.
rdat2 <- left_join(rdat, dat_exclusions, by = join_by(sub)) %>%
  filter(
    # Sufficient anat data
    anat_include == "anat_include",
    # Sufficient tfMRI or rsfMRI data
    rstask_include == "rstask_include" | rs_include == "rs_include",
  ) %>%
  select(-anat_include)

# Should be 1087

lefties <- rdat2 %>%
  filter(
    hand1 == "left"
  )

# Unrelated righties - randomly pick 100 for balanced sampling
righties <- rdat2 %>%
  filter(
    hand1 == "right",
    !(family %in% lefties$family)
  ) %>%
  slice_sample(n = nrow(lefties))

list1 <- bind_rows(righties, lefties) %>%
  select(sub, hand1) %>%
  arrange(sub)

leftovers <- rdat2 %>%
  filter(
    !(sub %in% list1$sub)
  ) %>%
  mutate(
    group = row_number() %/% 240
  ) %>%
  nest(.by = group)

write_csv(list1, here("list1.csv"))
write_csv(leftovers$data[[1]], here("list2.csv"))
write_csv(leftovers$data[[2]], here("list3.csv"))
write_csv(leftovers$data[[3]], here("list4.csv"))
write_csv(leftovers$data[[4]], here("list5.csv"))


## HI1 ====

hi1 <- filter(rdat2, hand1 == "right")
# table(hi1$rstask_include, hi1$rs_include)

## HI2a ====

rdat2 %>%
  group_by(hand1, rstask_include, rs_include) %>%
  summarize(
    n = n()
  )


## HI2b ====

rdat2 %>%
  group_by(hand2, rstask_include, rs_include) %>%
  summarize(
    n = n()
  )

rm(dat, dat_exclusions, dat0, rdat, rdat0, hi1)

# Do the balancing ====

# These are the overall ratios to check the folds and splits against

# Georgetown's ZIP code
set.seed(20057)

ratios <- list()
ratios["male"] <- sum(rdat2$gender == "M") / nrow(rdat2)
ratios["hand1"] <- sum(rdat2$hand1 == "right") / nrow(rdat2)
ratios["hand2"] <- sum(rdat2$hand2 == "right2") / nrow(rdat2)
mean_age <- mean(rdat2$age)

# Create the five mutually exclusive folds
folds <- rdat2 %>%
  select(-anat_include, -rs_include, -ehi) %>%
  mutate(
    fold = cluster_ra(clusters = family, conditions = LETTERS[1:5])
  ) %>%
  group_by(fold) %>%
  nest(.key = "test") %>%
  arrange(fold) %>%
  mutate(
    size = map_int(test, nrow)
  )

# Stop here and check if the sizes are unequally distributed
size_chisq <- chisq.test(tibble(folds$size, nrow(rdat2) / 5),
                         simulate.p.value = TRUE)

# This should be > .05
size_chisq$p.value

get_complement <- function(src, exclude_subs) {

  # This function gets all the rows that aren't included in the test sample

  result <- src %>%
    filter(
      !(sub %in% exclude_subs)
    )

  return(result)

}

folds_tt <- folds %>%
  select(-size) %>%
  mutate(
    train = map(test, ~get_complement(rdat2, .x$sub))
  ) %>%
  pivot_longer(c(test, train), names_to = "tt_group", values_to = "data") %>%
  mutate(
    n = map_int(data, nrow),
    age = map_dbl(data, ~mean(.x$age)),
    male_pct = map_dbl(data, ~sum(.x$gender == "M")) / n,
    hand1R_pct = map_dbl(data, ~sum(.x$hand1 == "right")) / n,
    hand2R_pct = map_dbl(data, ~sum(.x$hand2 == "right2")) / n
  )

# Check that the ages and ratios in each group don't differ from the mean

age_chisq <- tibble(age = folds_tt$age, avg_age = mean_age) %>%
  chisq.test(simulate.p.value = TRUE)

gender_chisq <- tibble(m = folds_tt$male_pct, avg_m = ratios$male) %>%
  chisq.test(simulate.p.value = TRUE)

hand1R_chisq <- tibble(m = folds_tt$hand1R_pct, avg_hand1R = ratios$hand1) %>%
  chisq.test(simulate.p.value = TRUE)

hand2R_chisq <- tibble(m = folds_tt$hand2R_pct, avg_hand2R = ratios$hand2) %>%
  chisq.test(simulate.p.value = TRUE)

sapply(list(age_chisq, gender_chisq, hand1R_chisq, hand2R_chisq),
       function(x) x$p.value)

folds_tt_final <- folds_tt %>%
  select(fold, tt_group, data) %>%
  unnest(data)

# lefties
(1 - ratios$hand1) * c(217, 897)
(1 - ratios$hand2) * c(217, 897)

# Get 30 training subs from Fold A for testing purposes
# a_train_rh <- folds_tt_final %>%
#   filter(
#     fold == "A",
#     tt_group == "train",
#     hand1 == "right"
#   ) %>%
#   slice_sample(n = 30) %>%
#   select(sub, fold, tt_group, hand1)
#
# write_csv(a_fold_rh, "/Volumes/thufir/HCPYA/lists/A_train_RH_n30.csv")

a_test_rh <- folds_tt_final %>%
  filter(
    fold == "A",
    tt_group == "test",
    hand1 == "right"
  ) %>%
  slice_sample(n = 30) %>%
  select(sub, fold, tt_group, hand1)

write_csv(a_test_rh, "/Volumes/thufir/HCPYA/lists/A_test_RH_n30.csv")
