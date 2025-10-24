library(tidyverse)
library(qs)
library(reticulate)

select <- dplyr::select

setwd("~/Google Drive/My Drive/Projects/hemisphere_fingerprinting/code/modeling/")

# demographic data ====

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
    EHI = Handedness
  ) %>%
  mutate(
    sub = paste0("sub-", Subject),
    handedness = if_else(EHI > 0, "righty", "lefty"),
    handedness2 = if_else(EHI >= 30, "righty2", "nonrighty")
  ) %>%
  select(sub, EHI, starts_with("handedness"))

ggplot(hands, aes(x = EHI)) +
  geom_histogram(bins = 40, boundary = 100) +
  scale_x_continuous(breaks = seq(-100, 100, 20)) +
  theme_bw()

# Imaging results =====

pconns <- list.files("../../data/pconns_0311_complete/", full.names = TRUE)

pconns_language <- list.files("../../data/pconns_0424_language/",
                              full.names = TRUE)

labels <- read_table("../../data/glasserhcp/labels.txt",
                     show_col_types = FALSE,
                     col_names = c("index", "label", "R", "G", "B", "A")) %>%
  select(index, label) %>%
  mutate(
    label = str_remove(label, "_ROI"),
    label2 = str_remove(label, "_")
  )

# This is all the combinations, to avoid 'X_Y' and 'Y_X' duplications
all_conns <- t(combn(labels$label2, 2)) %>%
  as_tibble() %>%
  mutate(
    conn = paste(V1, V2, sep = "_")
  )

## Function to read data ====

read_data <- function(path, labels, which) {

  pairwise <- read_table(path, col_names = labels, show_col_types = FALSE)

  # Remove upper triangle of duplicated r values
  pairwise[upper.tri(pairwise, diag = TRUE)] <- NA

  # Add ROI1 columns
  pairwise2 <- pairwise %>%
    add_column(ROI1 = labels, .before = 1)

  # Pivot longer and drop upper triangle
  pairwise3 <- pairwise2 %>%
    pivot_longer(-ROI1, names_to = "ROI2", values_to = "r") %>%
    na.omit()

  if (which == "hemi") {

    # Keep only within-hemisphere connections

    keep <- pairwise3 %>%
      filter(
        str_extract(ROI1, "^.") == str_extract(ROI2, "^.")
      ) %>%
      mutate(
        new_label = paste0(ROI1, "_", str_remove(ROI2, "^._"))
      ) %>%
      separate_wider_delim(new_label, delim = "_", names = c("hemi", "rois"),
                           too_many = "merge") %>%
      select(hemi, rois, r) %>%
      mutate(
        hemi = paste0(hemi, "H")
      ) %>%
      pivot_wider(names_from = rois, values_from = r)

  } else if (which == "trans") {

    # Keep only L<->R connections

    keep <- pairwise3 %>%
      filter(
        ROI1 != ROI2,
        str_extract(ROI1, "^.") != str_extract(ROI2, "^.")
      ) %>%
      mutate(
        new_label = paste0(str_remove(ROI1, "_"), "_", str_remove(ROI2, "_"))
      ) %>%
      select(-starts_with("ROI")) %>%
      pivot_wider(names_from = new_label, values_from = r)

  } else if (which == "all") {

    # Keep all unique values

    keep <- pairwise3 %>%
      filter(
      ) %>%
      mutate(
        new_label = paste0(ROI1, "_", ROI2)
      ) %>%
      select(-starts_with("ROI")) %>%
      pivot_wider(names_from = new_label, values_from = r)

  }

  return(keep)

}

subs_to_keep_rest <- unique(str_remove(basename(pconns), "_.*.txt"))
subs_to_keep_lang <- unique(str_remove(basename(pconns_language), "_.*.txt"))

## Create groups based on resting-state data ====

# These groups are carried forward into all manipulations

groups <- hands %>%
  distinct() %>%
  group_by(handedness) %>%
  mutate(
    group = cut_interval(row_number(), n = 5, labels = LETTERS[1:5])
  ) %>%
  ungroup() %>%
  select(sub, group, handedness)

groups_summary <- groups %>%
  group_by(group) %>%
  summarize(
    n_righty = sum(handedness == "righty"),
    n_lefty = sum(handedness == "lefty")
  )

# Hemi/trans connectome ====

info_cols <- left_join(hands, sex, by = join_by(sub)) %>%
  left_join(groups, by = join_by(sub, handedness))

## Hemiconnectome ====

if (!file.exists("inputs/hemiconnectome.rds")) {

  # Creates all hemiconnectomes, full size = length(p)*2
  all_data <- tibble(f = pconns) %>%
    mutate(
      pconn = map(f, ~read_data(.x, labels$label, which = "hemi"),
                  .progress = TRUE)
    )

  all_data2 <- all_data %>%
    unnest(pconn)

  all_data3 <- all_data2 %>%
    mutate(
      sub = basename(f) %>%
        str_remove("_task-rest.pconn.txt")
    ) %>%
    left_join(info_cols, by = join_by(sub)) %>%
    mutate(
      class = paste0(handedness, hemi),
      class2 = paste0(handedness2, hemi)
    ) %>%
    select(sub, gender, age_group, hemi, starts_with("handedness"),
           starts_with("class"), everything(), -f) %>%
    mutate(
      # Fisher Z transform correlation values - approximate normality
      across(c(where(is.numeric), -EHI), DescTools::FisherZ)
    )

  qsave(all_data3, "inputs/hemiconnectome.rds")
  py_save_object(all_data3, "inputs/hemiconnectome.pickle")

}

## Transconnectome ====

if (!file.exists("inputs/transconnectome.rds")) {

  # Creates all transconnectomes, full size = length(p)
  all_trans <- tibble(f = pconns) %>%
    mutate(
      pconn = map(f, ~read_data(.x, labels$label, which = "trans"),
                  .progress = TRUE)  )


  all_trans2 <- all_trans %>%
    unnest(pconn)

  all_trans3 <- all_trans2 %>%
    mutate(
      sub = basename(f) %>%
        str_remove("_task-rest.pconn.txt")
    ) %>%
    left_join(info_cols) %>%
    select(sub, gender, age_group, handedness, everything(), -f) %>%
    mutate(
      # Fisher Z transform correlation values - approximate normality
      across(c(where(is.numeric), -EHI), DescTools::FisherZ)
    )

  qsave(all_trans3, "inputs/transconnectome.rds")
  py_save_object(all_trans3, "inputs/transconnectome.pickle")

}

## Language data ====

  if (!file.exists("inputs/hemiconnectome_language.rds")) {

  all_data_lang <- tibble(f = pconns_language) %>%
    mutate(
      pconn = map(f, ~read_data(.x, labels$label, which = "hemi"),
                  .progress = TRUE)
    )

  all_trans_lang <- tibble(f = pconns_language) %>%
    mutate(
      pconn = map(f, ~read_data(.x, labels$label, which = "trans"),
                  .progress = TRUE)
    )

  all_lang_data2 <- all_data_lang %>%
    unnest(pconn)

  all_lang_data3 <- all_lang_data2 %>%
    mutate(
      sub = basename(f) %>%
        str_remove("_task-.*.pconn.txt")
    ) %>%
    left_join(hands, by = join_by(sub)) %>%
    left_join(sex, by = join_by(sub)) %>%
    mutate(
      class = paste0(handedness, hemi),
      class2 = paste0(handedness2, hemi)
    ) %>%
    select(sub, gender, age_group, hemi, starts_with("handedness"),
           starts_with("class"), everything(), -f) %>%
    left_join(groups, by = join_by(sub)) %>%
    mutate(
      # Fisher Z transform correlation values - approximate normality
      across(c(where(is.numeric), -EHI), DescTools::FisherZ)
    )

  qsave(all_lang_data3, "inputs/hemiconnectome_language.rds")
  py_save_object(all_lang_data3, "inputs/hemiconnectome_language.pickle")

}

if (!file.exists("inputs/transconnectome_language.rds")) {

  all_lang_trans2 <- all_trans_lang %>%
    unnest(pconn)

  all_lang_trans3 <- all_lang_trans2 %>%
    mutate(
      sub = basename(f) %>%
        str_remove("_task-rest.pconn.txt")
    ) %>%
    left_join(hands, by = join_by(sub)) %>%
    left_join(sex, by = join_by(sub)) %>%
    select(sub, gender, age_group, handedness, everything(), -f) %>%
    left_join(groups, by = join_by(sub)) %>%
    mutate(
      # Fisher Z transform correlation values - approximate normality
      across(c(where(is.numeric), -EHI), DescTools::FisherZ)
    )


  transconnectome_lang_Z <- all_lang_trans3 %>%
    mutate(
      # Fisher Z transform correlation values - approximate normality
      across(c(where(is.numeric), -EHI), DescTools::FisherZ)
    )

  qsave(transconnectome_lang_Z, "inputs/transconnectome_language.rds")
  py_save_object(transconnectome_lang_Z,
                 "inputs/transconnectome_language.pickle")

}

## Full connectome ====

if (!file.exists("inputs/connectome.rds")) {

  # Creates all complete connectomes, full size = length(p)
  all_complete <- tibble(f = pconns) %>%
    mutate(
      pconn = map(f, ~read_data(.x, labels$label2, which = "all"),
                  .progress = TRUE)
    )

  all_complete2 <- all_complete %>%
    unnest(pconn)

  all_complete3 <- all_complete2 %>%
    mutate(
      sub = basename(f) %>%
        str_remove("_task-rest.pconn.txt")
    ) %>%
    left_join(hands, by = join_by(sub)) %>%
    left_join(sex, by = join_by(sub)) %>%
    select(sub, gender, age_group, handedness, everything(), -f) %>%
    left_join(groups, by = join_by(sub)) %>%
    mutate(
      # Fisher Z transform correlation values - approximate normality
      across(c(where(is.numeric), -EHI), DescTools::FisherZ)
    )

  qsave(all_complete3, "inputs/connectome.rds")
  py_save_object(all_complete3, "inputs/connectome.pickle")

}

## Hemireversome ====

if (!file.exists("inputs/reversed.rds")) {

  reversed_pconns <- list.files("../../data/flipped_pconns/",
                                pattern = "*.pconn.txt",
                                full.names = TRUE)

  # Creates all complete connectomes, full size = length(p)
  rev <- tibble(f = reversed_pconns) %>%
    mutate(
      pconn = map(f, ~read_data(.x, labels$label, which = "hemi"),
                  .progress = TRUE)
    )

  rev2 <- rev %>%
    mutate(
      sub = str_extract(f, "sub-[0-9]*")
    ) %>%
    select(sub, pconn) %>%
    unnest(pconn)

  rev_final <- rev2 %>%
    left_join(info_cols, join_by(sub)) %>%
    select(sub, gender, age_group, EHI, group, starts_with("handedness"),
           everything()) %>%
    mutate(
      # Fisher Z transform correlation values - approximate normality
      across(c(where(is.numeric), -EHI), DescTools::FisherZ)
    )

  qsave(rev_final, "inputs/reversed.rds")
  py_save_object(rev_final, "inputs/reversed.pickle")

}

