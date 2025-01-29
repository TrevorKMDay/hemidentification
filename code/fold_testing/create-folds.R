library(tidyverse)

setwd(paste0("~/Library/CloudStorage/GoogleDrive-td758@georgetown.edu/",
             "My Drive/Projects/hemisphere_fingerprinting/"))

# Load data ====

dat0 <- read_csv("data/unrestricted_data.csv", show_col_types = FALSE)  %>%
  select_all(~str_remove_all(., "-|3T_"))

rdat0 <- read_csv("data/restricted_data.csv", show_col_types = FALSE)

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
  ) %>%
  group_by(family) %>%
  mutate(
    fam_n = n(),
    fam_m_age = mean(age)
  )

rdat2 <- left_join(rdat, dat_exclusions, by = join_by(sub)) %>%
  filter(
    # Sufficient anat data
    anat_include == "anat_include",
    # Sufficient tfMRI or rsfMRI data
    rstask_include == "rstask_include" | rs_include == "rs_include",
  )

# Should be 1087

## HI1 ====

hi1 <- filter(rdat2, hand1 == "right")
table(hi1$rstask_include, hi1$rs_include)

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

# Do the balancing

rdat2_fams <- rdat2 %>%
  group_by(family, fam_n, fam_m_age) %>%
  nest() 

create_split <- function(fams) {
  
  n1 <- nrow(fams)
  ntrain <- round(0.8 * n1)
  
  row <- 1:n1
  training <- sample.int(n1, ntrain)
  test <- row[!(row %in% training)]
  
  train_set <- unnest(fams[training, ], data)
  test_set <- unnest(fams[test, ], data)
  
  return(tibble(train = list(train_set), test = list(test_set)))
  
}

random_splits <- tibble(f = 1:100) %>%
  rowwise() %>%
  mutate(
    split = list(create_split(rdat2_fams))
  ) %>%
  unnest(split) %>%
  mutate(
    train_ratio = map2_dbl(train, test, ~ nrow(.x) / (nrow(.x) + nrow(.y))),
    
    train_men = map_dbl(train, ~sum(.x$gender == "M") / nrow(.x)),
    test_men = map_dbl(test, ~sum(.x$gender == "M") / nrow(.x)),
    
    train_mean_age = map_dbl(train, ~mean(.x$age)),
    test_mean_age = map_dbl(test, ~mean(.x$age)),
    
    train_rhand1 = map_dbl(train, ~sum(.x$hand1 == "right") / nrow(.x)),
    test_rhand1 = map_dbl(test, ~sum(.x$hand1 == "right") / nrow(.x)),
    
    train_rhand2 = map_dbl(train, ~sum(.x$hand2 == "right2") / nrow(.x)),
    test_rhand2 = map_dbl(test, ~sum(.x$hand2 == "right2") / nrow(.x)),
    
  )


rdat2_male_ratio <- sum(rdat2$gender == "M") / nrow(rdat2)
rdat2_hand1_ratio <- sum(rdat2$hand1 == "right") / nrow(rdat2)
rdat2_hand2_ratio <- sum(rdat2$hand2 == "right2") / nrow(rdat2)

balanced_splits <- random_splits %>%
  mutate(
    age_error = abs(test_mean_age - train_mean_age) / ((test_mean_age + train_mean_age) / 2)
  ) %>%
  filter(
    
    abs(train_ratio - 0.8) < .02,
    age_error < .02,
    
    abs(train_men - rdat2_male_ratio) < .02,
    abs(test_men  - rdat2_male_ratio) < .02,
    
    abs(train_rhand1 - rdat2_hand1_ratio) < .02,
    abs(test_rhand1 - rdat2_hand1_ratio) < .02,
    
    abs(train_rhand2 - rdat2_hand2_ratio) < .02,
    abs(test_rhand2 - rdat2_hand2_ratio) < .02
    
  )

