library(tidyverse)
library(qs2)
library(reticulate)

select <- dplyr::select

setwd(str_glue("~/Google Drive/My Drive/Projects/hemisphere_fingerprinting/",
               "code/modeling_new/"))

# Demographic data ====

info_cols <- read_csv("folds_to_use2.csv", show_col_types = FALSE) 

# Imaging results =====

pconns <- list.files("../../data/pconns_0311_complete/", full.names = TRUE)

pconns_language <- list.files("../../data/pconns_0424_language/",
                              full.names = TRUE)

pconns_LL <- list.files("../../data/glasserLL/", full.names = TRUE)
pconns_RR <- list.files("../../data/glasserRR/", full.names = TRUE)

labels <- read_table("../../data/glasserhcp/labels.txt",
                     show_col_types = FALSE,
                     col_names = c("index", "label", "R", "G", "B", "A")) %>%
  select(index, label) %>%
  mutate(
    label = str_remove(label, "_ROI"),
    label2 = str_remove(label, "_")
  )

flipped_labels <- read_table("flipped_labels_0to760.txt",
                     show_col_types = FALSE,
                     col_names = c("index", "label", "R", "G", "B", "A")) %>%
  select(index, label) %>%
  filter_out(
    label %in% c("INDEXMAX", "???")
  ) %>%
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

# Networks ====

networks <- read_csv("../../data/ColeAnticevicNetPartition/all_parcels_to_networks.csv",
                     show_col_types = FALSE)

motor_rois <- c(networks$label2[networks$L == "Somatomotor"],
                networks$label2[networks$R == "Somatomotor"]) %>%
  unique()

# Function to read data ====

read_data <- function(path, columns, which, filter_missing_rois = FALSE) {

  pairwise <- read_table(path, col_names = columns, show_col_types = FALSE,
                          na = "nan")
  
  # For the matrices where ROIs are missing, filter rows/cols that are all
  # missing values, except the `1` on the diagonals
  if (filter_missing_rois)
    pairwise <- pairwise[rowSums(!is.na(pairwise)) > 1, 
                         colSums(!is.na(pairwise)) > 1]

  probs <- problems(pairwise)

  if (nrow(probs) > 0) {
    print("Parsing failure")
    return(NA)
  } else {

    # Remove upper triangle of duplicated r values
    pairwise[upper.tri(pairwise, diag = TRUE)] <- NA

    # Add ROI1 columns
    pairwise2 <- pairwise %>%
      add_column(ROI1 = colnames(pairwise), .before = 1)

    # Pivot longer and drop upper triangle
    pairwise3 <- pairwise2 %>%
      pivot_longer(-ROI1, names_to = "ROI2", values_to = "r") %>%
      na.omit()

    if (which == "hemi") {

      # Keep only within-hemisphere connections

      keep <- pairwise3 %>%
        filter(
          str_extract(ROI1, "^.*_") == str_extract(ROI2, "^.*_")
        ) %>%
        mutate(
          new_label = paste0(ROI1, "_", str_remove(ROI2, "^.*_"))
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

}

subs_to_keep_rest <- unique(str_remove(basename(pconns), "_.*.txt"))
subs_to_keep_lang <- unique(str_remove(basename(pconns_language), "_.*.txt"))

# Create datasets ====

## Hemiconnectome ====

# Use parallelism to read in data, otherwise it's quite long
library(furrr)
future::plan(strategy = multicore, workers = 11)

if (!file.exists("inputs/hemiconnectome.rds")) {

  # Creates all hemiconnectomes, full size = length(p)*2
  all_data <- tibble(f = pconns) %>%
    mutate(
      pconn = future_map(f, ~read_data(.x, labels$label, which = "hemi"),
                         .progress = TRUE)
    )
  
  # Add cx-prefix to easily identify columns
  all_data2 <- all_data %>%
    unnest(pconn) %>%
    rename_at(vars(contains("_")), ~str_c("cx_", .)) 

  all_data3 <- all_data2 %>%
    mutate(
      sub = basename(f) %>%
        str_remove("_task-rest.pconn.txt"),
    ) %>%
    filter(
      sub %in% info_cols$sub
    ) %>%
    inner_join(info_cols, ., by = join_by(sub)) %>%
    mutate(
      class = paste0(handedness, hemi),
    ) %>%
    select(
      # Reorder cols and drop exact age
      sub, group, hemi_fold, gender, age_group, hemi,
      starts_with("handedness"), starts_with("class"),
      everything(),
      -c(f, age)
    ) %>%
    mutate(
      # Fisher Z transform correlation values - approximate normality
      across(c(where(is.numeric), -EHI), DescTools::FisherZ)
    )

  qs_save(all_data3, "inputs/hemiconnectome.rds", nthreads = 10)
  py_save_object(all_data3, "inputs/hemiconnectome.pickle")
  
}

## Normalized hemiconnectome ====

if (!file.exists("inputs/hc_hpZ.rds")) {

  hc <- qs_read("inputs/hemiconnectome.rds", nthreads = 10)

  hc_long <- hc %>%
    pivot_longer(starts_with("cx_"))

  # Normalize within-hemisphere
  hc_standardized <- hc_long %>%
    group_by(sub, hemi) %>%
    mutate(
      value = scale(value)[, 1]
    )
 
  hc_hpZ <- pivot_wider(hc_standardized)
  
  qs_save(all_data3, "inputs/hc_hpZ.rds", nthreads = 10)
  py_save_object(all_data3, "inputs/hc_hpZ.pickle")

}

## Language data ====

if (!file.exists("inputs/hemiconnectome_language_hpZ.rds")) {

  all_data_lang <- tibble(f = pconns_language) %>%
    mutate(
      pconn = map(f, ~read_data(.x, labels$label, which = "hemi"),
                  .progress = TRUE)
    )

  all_lang_data2 <- all_data_lang %>%
    unnest(pconn) %>%
    rename_at(vars(contains("_")), ~str_c("cx_", .)) 

  all_lang_data3 <- all_lang_data2 %>%
    mutate(
      sub = basename(f) %>%
        str_remove("_task-.*.pconn.txt")
    ) %>%
    inner_join(info_cols, ., by = join_by(sub)) %>%
    mutate(
      class = paste0(handedness, hemi),
    ) %>%
    select(sub, gender, age_group, hemi, starts_with("handedness"),
           starts_with("class"), everything(), -f, -age) %>%
    mutate(
      # Fisher Z transform correlation values - approximate normality
      across(c(where(is.numeric), -EHI), DescTools::FisherZ)
    )

  lang_long <- all_lang_data3 %>%
    pivot_longer(starts_with("cx_"))

  # Normalize within-hemisphere
  lang_standardized <- lang_long %>%
    group_by(sub, hemi) %>%
    mutate(
      value = scale(value)[, 1]
    )
 
  xx <- lang_standardized %>%
    ungroup() %>%
    dplyr::summarise(
      n = dplyr::n(), 
      .by = c(sub, gender, age_group, hemi, handedness, class, EHI, hemi_fold, group,
    name)) %>%
    dplyr::filter(n > 1L) 
  
  lang_hpZ <- pivot_wider(ungroup(lang_standardized))

  qs_save(all_lang_data3, "inputs/hemiconnectome_language_hpZ.rds")
  py_save_object(all_lang_data3, 
                 "inputs/hemiconnectome_language_hpZ.pickle")

}

## LL/RR ====

flips <- tribble(
    ~name, ~pconns,
    "LL", pconns_LL,
    "RR", pconns_RR
  )

for (i in 1:nrow(flips)) {

  label <- flips$name[i]
  pconns <- flips$pconns[[i]]

  print(label)

  if (!file.exists(str_glue("inputs/hemiconnectome_{label}.rds"))) {

    # Creates all hemiconnectomes, full size = length(p)*2
    all_data <- tibble(f = pconns) %>%
      mutate(
        pconn = map(f, ~read_data(.x, flipped_labels$label, 
                                         which = "hemi", 
                                         filter_missing_rois = TRUE),
                          .progress = TRUE)
      )
    
    # Add cx-prefix to easily identify columns
    all_data2 <- all_data %>%
      unnest(pconn) %>%
      rename_at(vars(contains("_")), ~str_c("cx_", .)) 

    all_data3 <- all_data2 %>%
      mutate(
        sub = basename(f) %>%
          str_remove("_task-rest.*.pconn.txt"),
      ) %>%
      filter(
        sub %in% info_cols$sub
      ) %>%
      inner_join(info_cols, ., by = join_by(sub)) %>%
      mutate(
        class = paste0(handedness, hemi),
      ) %>%
      select(
        # Reorder cols and drop exact age
        sub, group, hemi_fold, gender, age_group, hemi,
        starts_with("handedness"), starts_with("class"),
        everything(),
        -c(f, age)
      ) %>%
      mutate(
        # Fisher Z transform correlation values - approximate normality
        across(c(where(is.numeric), -EHI), DescTools::FisherZ)
      )
    
    hc_long <- all_data3 %>%
      pivot_longer(starts_with("cx_"))

    # Normalize within-hemisphere
    hc_standardized <- hc_long %>%
      group_by(sub, hemi) %>%
      mutate(
        value = scale(value)[, 1]
      )
    
    hc_hpZ <- pivot_wider(hc_standardized)

    qs_save(hc_hpZ, str_glue("inputs/hemiconnectome_{label}.rds"), 
            nthreads = 10)
    
    py_save_object(hc_hpZ, str_glue("inputs/hemiconnectome_{label}.pickle"))
    
  }
  
}

## Transconnectome ====

if (!file.exists("inputs/transconnectome.rds")) {

  # Creates all transconnectomes, full size = length(p)
  all_trans <- tibble(f = pconns) %>%
    mutate(
      pconn = map(f, ~read_data(.x, labels$label, which = "trans"),
                  .progress = TRUE))

  all_trans2 <- all_trans %>%
    unnest(pconn) %>%
    rename_at(vars(contains("_")), ~str_c("cx_", .)) 

  all_trans3 <- all_trans2 %>%
    mutate(
      sub = basename(f) %>%
        str_remove("_task-rest.pconn.txt")
    ) %>%
    left_join(info_cols, by = join_by(sub)) %>%
    select(sub, gender, age_group, handedness, everything(), -f, -age) %>%
    mutate(
      # Fisher Z transform correlation values - approximate normality
      across(c(where(is.numeric), -EHI), ~scale(DescTools::FisherZ(.x))[, 1])
    )
  
  # check for missing vlaues
  which(colSums(is.na(all_trans3)) > 0)

  qs_save(all_trans3, "inputs/transconnectome.rds")
  py_save_object(all_trans3, "inputs/transconnectome.pickle")

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
    unnest(pconn)  %>%
    rename_at(vars(contains("_")), ~str_c("cx_", .)) 

  all_complete3 <- all_complete2 %>%
    mutate(
      sub = basename(f) %>%
        str_remove("_task-rest.pconn.txt")
    ) %>%
    left_join(info_cols) %>%
    select(sub, gender, age_group, handedness, everything(), -f, -age) %>%
    mutate(
      # Fisher Z transform correlation values - approximate normality
      # then scale columns
      across(c(where(is.numeric), -EHI), ~scale(DescTools::FisherZ(.x))[, 1])
    )

  which(colSums(is.na(all_complete3)) > 0)
  
  qs_save(all_complete3, "inputs/connectome.rds")
  py_save_object(all_complete3, "inputs/connectome.pickle")

}

## Oversample HC =====

library(groupdata2)

if (!file.exists("hc_oversampled.rds")) {

  hc <- qs_read("inputs/hc_hpZ.rds")

  hc_oversampled <- upsample(hc, cat_col = "class")

  qs_save(hc_oversampled, "inputs/hc_oversampled.rds")
  py_save_object(hc_oversampled, "inputs/hc_oversampled.pickle")

}