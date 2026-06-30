library(tidyverse)

# Parallel setup ====

library(parallel)
library(doParallel)
library(foreach)

cores <- detectCores()
cl <- makeCluster(cores[1] - 1) #not to overload your computer
registerDoParallel(cl)

# Setup

select <- dplyr::select
setwd("~/MyDrive/Projects/hemisphere_fingerprinting/code/modeling_new")

bs_reps <- 1000

# Code to bootstrap MCC 

source("../bootstrap_mcc.R")

# Read in results

folds <- read_csv("folds_to_use2.csv", show_col_types = FALSE) %>%
  select(sub, group)

result_files <- Sys.glob("03-*/results/*/*.csv")

results0 <- tibble(f = result_files) %>%
  mutate(
    bn = basename(f),
    data = map(f, read_csv, show_col_types = FALSE),
    hand = str_extract(f, "[a-z]*_handed")
  ) %>%
  select(-f) %>%
  separate_wider_delim(bn, delim = "_", 
                       names = c("input", "clf", "outcome", NA)) %>%
  unnest(data) %>%
  select(-`...1`)

overall_file <- str_glue("mccresults_reps-{bs_reps}.rds")

if (!file.exists(overall_file)) { 
  
  # MCC overall 
  results_overall0 <- results0 %>%
    group_by(hand, input, clf, outcome) %>%
    nest() %>%
    mutate(

      # Parallel bootstrapping done here:
      mcc = map(data, ~bootstrap_mcc(.x, R = bs_reps, ground = "hemi"), 
                .progress = TRUE)

    )
  
  qs2::qs_save(results_overall0, overall_file)

} else {

  results_overall0 <- qs2::qs_read(overall_file)

}

results_overall <- results_overall0 %>%
  unnest(mcc) %>%
  mutate(
    mcc_ci_hi = if_else(mcc_ci_hi > 1, 1, mcc_ci_hi),
    input = str_remove(input, "data-")
  )

bygroup_file <- str_glue("mccresults_reps-{bs_reps}_bygroup.rds")

if (!file.exists(bygroup_file)) { 

  # Get MCC by group
  results_bygroup0 <- results0 %>%
    left_join(folds, by = join_by(sub)) %>%
    group_by(hand, input, clf, outcome, group) %>%
    nest() %>%
    mutate(

      mcc = map(data, ~bootstrap_mcc(.x, R = bs_reps, ground = "hemi"),
                .progress = TRUE)

    )
  
  qs2::qs_save(results_bygroup0, bygroup_file)

} else {

  results_bygroup0 <- qs2::qs_read(bygroup_file)
  
}

results_bygroup <- results_bygroup0 %>%
  unnest(mcc) %>%
  mutate(
    mcc_ci_hi = if_else(mcc_ci_hi > 1, 1, mcc_ci_hi),
    input = str_remove(input, "data-")
  )

# Put together for plotting 

lda_results_overall <-  results_overall %>%
  filter(
    clf == "clf-lda",
  ) %>%
  filter_out(
    input %in% c("hc.RHD", "hemiconnectome.LHD") |
    str_detect(input, "^hc.1p.")
  ) %>%
  mutate(

    input_data = str_extract(input, "(hc|lang).hpZ.") %>%
      replace_values(
        "hc.hpZ." ~ "HC",
        "lang.hpZ." ~ "Lang" 
      ),

    one_per = if_else(str_detect(input, "1p"), " (1 per)", ""),

    hand_tag = str_extract(input, "[LR]HD.*") %>%
      str_replace("rdm", " (random)") %>%
      str_replace("prd", " (paired)") %>%
      str_replace("LHD", "Sinistral") %>%
      str_replace("RHD", "Dextral"),

    name = str_glue("{input_data}{one_per}")

  )


lda_results_bygroup <- results_bygroup %>%
  filter(
    clf == "clf-lda",
  ) %>%
  filter_out(
    input %in% c("hc.RHD", "hemiconnectome.LHD") |
    str_detect(input, "^hc.1p.")
  ) %>%
  mutate(

    input_data = str_extract(input, "(hc|lang).hpZ.") %>%
      replace_values(
        "hc.hpZ." ~ "HC",
        "lang.hpZ." ~ "Lang" 
      ),

    one_per = if_else(str_detect(input, "1p"), " (1 per)", ""),

    hand_tag = str_extract(input, "[LR]HD.*") %>%
      str_replace("rdm", " (random)") %>%
      str_replace("prd", " (paired)") %>%
      str_replace("LHD", "Sinistral") %>%
      str_replace("RHD", "Dextral"),

    name = str_glue("{input_data}{one_per}")

  )

ggplot(lda_results_bygroup, aes(x = hand_tag, color = group)) +
  geom_hline(yintercept = 0, color = "red") +
  geom_pointrange(aes(y = mcc_mean, ymin = mcc_ci_lo, ymax = mcc_ci_hi),
                  position = position_dodge2(width = .25),
                  alpha  = .5) +
  geom_pointrange(
    data = lda_results_overall,
    aes(y = mcc_mean, ymin = mcc_ci_lo, ymax = mcc_ci_hi),
    color = "black"
  ) +
  scale_x_discrete(
    limits = c("Dextral", "Sinistral", "Dextral (random)", "Dextral (paired)"),
    labels = c("Dextral", "Sinistral", "Dextral\n(random)", "Dextral\n(paired)")
  ) +
  scale_y_continuous(limits = c(NA, 1)) +
  facet_wrap(vars(name)) +
  theme_bw() +
  labs(x = "Data", y = "MCC (95% CI)", color = "Fold") 

ggsave("plots/mcc_results_singlehandedness.png", width = 6.5, height = 6)

# Run test

lda_results_test <- lda_results_bygroup %>%
  ungroup() %>%
  select(input, input_data, one_per, group, hand_tag, mcc_mean) %>%
  mutate(

    one_per = one_per != "",

    hand_tag = replace_values(
      hand_tag,
      "RHD (paired)" ~ "RHDpaired",
      "RHD (random)" ~ "RHDrandom"
    )

  )

mcc_model <- aov(mcc_mean ~ group + input_data + one_per + hand_tag, 
                data = lda_results_test) 

summary(mcc_model)

mcc_model_tukey <- TukeyHSD(mcc_model, which = "hand_tag")

mcc_model_tukey_printable <- as_tibble(mcc_model_tukey$hand_tag, 
                                       rownames = "pair") %>%
  mutate(
    beta = str_glue("{round(diff, 3)} [{round(lwr, 3)}, {round(upr, 3)}]"),
    p_adj = round(`p adj`, 3)
  ) %>%
  select(pair, beta, p_adj)

clipr::write_clip(mcc_model_tukey_printable)
