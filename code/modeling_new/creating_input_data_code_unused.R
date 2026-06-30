

## Language data ====

  

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

  qs_save(transconnectome_lang_Z, "inputs/transconnectome_language.rds")
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

  all_complete3 %>%
    select(sub, matches("^[RL]")) %>%
    pivot_longer(-sub) %>%
    group_by(sub) %>%
    summarize(
      mean = mean(value),
      sd = sd(value)
    )


  qs_save(all_complete3, "inputs/connectome.rds")
  py_save_object(all_complete3, "inputs/connectome.pickle")

  # Motor

  motor_regex <- paste0(motor_rois, collapse = "|")

  all_complete3 <- qs::qread("inputs/connectome.rds", nthreads = 10)

  motor <- all_complete3 %>%
    select(sub, gender, group, age_group, handedness,
            matches(str_glue("^[RL]({motor_regex})_[RL]({motor_regex})$")))

  qs_save(motor, "inputs/motor_connectome.rds")
  py_save_object(motor, "inputs/motor_connectome.pickle")

  ## Z score =====

  ### All data =====

  all_complete_Z <- all_complete3 %>%
    pivot_longer(-c(sub, group, gender, age_group, starts_with("handedness"),
                    EHI),
                 values_to = "r")

  filter(all_complete_Z, abs(r) == 100)

  all_complete_Z2 <- all_complete_Z %>%
    group_by(sub) %>%
    mutate(
      rZ = scale(r)[, 1]
    ) %>%
    select(-r)

  all_complete_Z2 %>%
    select(sub, matches("^[RL]")) %>%
    pivot_longer(-sub) %>%
    group_by(sub) %>%
    summarize(
      mean = mean(value),
      sd = sd(value)
    )

  all_complete_Z3 <- all_complete_Z2 %>%
    pivot_wider(values_from = rZ, names_from = name)

  qs_save(all_complete_Z3, "inputs/connectome_Z.rds")
  py_save_object(all_complete_Z3, "inputs/connectome_Z.pickle")

  ### by hemisphere ====

  all_Z_byhemi <- all_complete_Z2 %>%
    mutate(
      h1 = str_extract(name, "^."),
      h2 = str_remove(str_extract(name, "_."), "_")
    ) %>%
    filter(
      h1 == h2
    ) %>%
    mutate(
      hemi = paste0(h1, "H"),
      name = str_remove(name, "^[LR]") %>%
        str_replace("_[LR]", "_")
    ) %>%
    select(-handedness2, -handedness.y, -h1, -h2)

  all_Z_byhemi %>%
    group_by(sub) %>%
    summarize(
      mean = mean(rZ),
      sd = sd(rZ)
    )

  all_Z_byhemi %>%
    group_by(sub,hemi) %>%
    summarize(
      mean = mean(rZ),
      sd = sd(rZ)
    )

  all_Z_byhemi_wide <- all_Z_byhemi %>%
    pivot_wider(names_from = name, values_from = rZ) %>%
    rename(
      handedness = handedness.x
    ) %>%
    mutate(
      class = paste0(handedness, hemi)
    )

  qs_save(all_Z_byhemi_wide, "inputs/hemiconnectome_Z.rds")
  py_save_object(all_Z_byhemi_wide, "inputs/hemiconnectome_Z.pickle")

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

  qs_save(rev_final, "inputs/reversed.rds")
  py_save_object(rev_final, "inputs/reversed.pickle")

}

# First half ====

if (!file.exists("inputs/hemiconnectome_run-1.rds")) {

  run1_pconns <- list.files("../../data/pconns_run-1_1113/",
                            pattern = ".*_seg-Glasser_.*.pconn.txt",
                            full.names = TRUE)

  # Creates all complete connectomes, full size = length(p)
  run1 <- tibble(f = run1_pconns) %>%
    mutate(
      pconn = map(f, ~read_data(.x, labels$label, which = "hemi"),
                  .progress = TRUE)
    )

  run1_2 <- run1 %>%
    mutate(
      sub = str_extract(f, "sub-[0-9]*")
    ) %>%
    select(sub, pconn) %>%
    unnest(pconn)

  run1_final <- run1_2 %>%
    left_join(info_cols, join_by(sub)) %>%
    select(sub, gender, age_group, EHI, group, starts_with("handedness"),
           everything()) %>%
    mutate(
      # Fisher Z transform correlation values - approximate normality
      across(c(where(is.numeric), -EHI), DescTools::FisherZ),
      class = paste0(handedness, hemi),
    ) %>%
    na.omit()

  # not sure why there's NAs - toss em

  qs_save(run1_final, "inputs/half_hemiconnectome.rds")
  py_save_object(run1_final, "inputs/half_hemiconnectome.pickle")

  # z score

  # run1_final_Z <- run1_final %>%
  #   pivot_longer(-c(hemi, class, colnames(info_cols))) %>%
  #   grou

}

# Randomly rearrange rows ======

# Create a consistent look up table to rename all ROIs to the same label in 
# each shuffle.

# Number of connections
n_conn <- (180 ** 2 - 180) / 2

set.seed(20007)

new_orders <- map(1:10, ~sample(1:n_conn, n_conn))

for (i in 1:10) {

  temp <- all_data3 %>%
    select(-class) %>%
    pivot_longer(-any_of(c(colnames(info_cols), "hemi"))) %>%
    pivot_wider(names_from = hemi) %>%
    group_by(sub) %>%
    mutate(
      # n = row_number()
    )

  reordered <- temp %>%
    mutate(

      # Sanity check this works right by shuffling n and verify each label
      # receives the same n
      # n = n[new_orders[[i]]]
      LH = LH[new_orders[[i]]],

      # Center each individual's hemiconnectome separately
      across(c(LH, RH), ~scale(.x)[,1])

    )
  
  reor_wide <- reordered %>%
    pivot_longer(c("RH", "LH"), names_to = "hemi") %>%
    pivot_wider()

  name <- str_glue("shuffle{str_pad(i, 2, 'left', '0')}")
  py_save_object(reor_wide, str_glue("inputs/shuffled_labels/{name}.pickle"))

  message(str_glue("Done with {name}!"))
 
}

# RH > LH
t.test(reordered$RH, reordered$LH, paired = TRUE)

# Shuffle independently

set.seed(20007)
shuffle2 <- all_data3 %>%
  select(-class) %>%
  pivot_longer(-any_of(c(colnames(info_cols), "hemi"))) %>%
  group_by(sub, hemi) %>%
  mutate(
    # Randomly reorder each individuals hemispheres separately
    value = sample(value, length(value))
  )

shuffle3 <- pivot_wider(shuffle2)
py_save_object(shuffle3, str_glue("inputs/shuffled_labels/shuffled_separately.pickle"))


# Completely random =====

set.seed(20007)
completely_random0 <- all_data3 %>%
  select(-class) %>%
  pivot_longer(-any_of(c(colnames(info_cols), "hemi"))) %>%
  group_by(sub) %>%
  mutate(
    value = sample(value, size = length(value))
  )

completely_random <- pivot_wider(completely_random0)
py_save_object(completely_random, str_glue("inputs/completely_random.pickle"))

# Summary values 

summary_values0 <- all_data3 %>%
  select(-class) %>%
  pivot_longer(-any_of(c(colnames(info_cols), "hemi"))) %>%
  group_by(group, sub, hemi)

summary_values <- summary_values0 %>%
  summarize(

    mean_edge_str = mean(value),
    median_edge_str = median(value),
    sumsq_edge_str = sum(value**2),
    var_edge_str = var(value)

  )

py_save_object(summary_values, str_glue("inputs/summary_values.pickle"))

# mean values 

mean_r <-  all_data3 %>%
  select(-class) %>%
  pivot_longer(-any_of(c(colnames(info_cols), "hemi"))) %>%
  group_by(name, hemi) %>%
  summarize(
    r = mean(value)
  ) %>%
  pivot_wider(names_from = hemi, values_from = r) %>%
  mutate(
    diff = LH - RH
  ) %>%
  arrange(desc(abs(diff)))
