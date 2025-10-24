library(tidyverse)
library(reticulate)
library(qs)
library(infer)
library(igraph)
library(patchwork)
library(Matrix)

# Setup ====

setwd("~/MyDrive/Projects/hemisphere_fingerprinting/code/enrichment/")
select <- dplyr::select

networks0 <- read_csv(file = paste0("../../data/ColeAnticevicNetPartition/",
                                    "all_parcels_to_networks.csv"),
                      show_col_types = FALSE)

nets <- networks0 %>%
  mutate(
    # Use LH network assignment, and combine Visual networks 1 and 2
    network = if_else(L == "Visual2", "Visual", L),
    net_short = case_match(
      network,
      "Auditory" ~ "AudN",
      "CinguloOpercular" ~ "CON",
      "Default" ~ "DMN",
      "DorsalAttention" ~ "DAN",
      "Frontoparietal" ~ "FPN",
      "Language" ~ "LN",
      "OrbitoAffective" ~ "OAN",
      "PosteriorMultimodal" ~ "pMN",
      "Somatomotor" ~ "SMN",
      "VentralMultimodal" ~ "vMN",
      "Visual" ~ "VisN"
    )
  ) %>%
  rename(
    roi = label2
  ) %>%
  select(roi, network, net_short)

nets_ordered <- nets %>%
  arrange(net_short, roi) %>%
  mutate(
    order = row_number(),
  ) %>%
  group_by(network) %>%
  mutate(
    order_in_net = row_number()
  )

write_rds(nets_ordered, "nets_ordered.rds")

# For plotting
nets_squares <- nets_ordered %>%
  filter(
    order_in_net == 1 | order_in_net == max(order_in_net)
  ) %>%
  mutate(
    startend = if_else(order_in_net == 1, "start", "end")
  ) %>%
  pivot_wider(id_cols = c(network, net_short), names_from = startend,
              values_from = roi) %>%
  mutate(
    # Replace 1 with 0 so it doesn't overlap the first ROI
    start_order = match(start, nets_ordered$roi),
    end_order = match(end, nets_ordered$roi)
  )

write_rds(nets_squares, "nets_squares.rds")

# Results ====

true_model <- read_csv("../modeling/lda_full/full_scalings.csv",
                       show_col_types = FALSE)

true_model_top12 <- true_model %>%
  arrange(desc(LD1)) %>%
  mutate(
    n = row_number()
  ) %>%
  filter(
    n <= 12
  )

second_model <- read_csv("../modeling/lda_full/full_class_scalings.csv",
                         show_col_types = FALSE)

cor(true_model$LD1, second_model$LD1) # 0.86
cor(true_model$LD1, second_model$LD1, method = "s") # 0.85

scnd_model_top12 <- second_model %>%
  arrange(desc(LD1)) %>%
  mutate(
    n = row_number()
  ) %>%
  filter(
    n <= 12
  )

intersect(true_model_top12$feature, scnd_model_top12$feature) # 6 / 12

## Overlap ====

intersect_features <- function(model1, model2, n) {

  model1 <- arrange(model1, desc(LD1))
  model2 <- arrange(model2, desc(LD1))

  overlap <- intersect(model1$feature[1:n], model2$feature[1:n])

  return(length(overlap))

}

n_overlapping <- tibble(n = seq(10, nrow(true_model), by = 200)) %>%
  mutate(
    intersect = map_int(n, ~intersect_features(true_model, second_model, .x)),
    pct = intersect / n
  )

ggplot(n_overlapping, aes(x = n, y = pct)) +
  geom_point() +
  theme_bw() +
  labs(title = "# overlapping features in top n features")

## Read in bootstrapping ====

bs_file <- "lda_full_bootstrapped.rds"
bs_file2 <- "lda_full_bootstrapped2.rds"
bs_fileL <- "lda_full_lefies.rds"

# See script 01a-bootstrap.R

bs <- qread(bs_file, nthreads = 2)
bs2 <- qread(bs_file2, nthreads = 10)
bsL <- qread(bs_fileL, nthreads = 10)

# Run tests ====

## 1-class =====

test_1class <- bs %>%
  group_by(feature) %>%
  nest() %>%
  left_join(true_model, join_by(feature)) %>%
  mutate(

    p = map2_dbl(data, LD1,
                  ~get_p_value(tibble(LD1 = .x$LD1), .y,
                               direction = "two_sided")[[1]],
                 .progress = TRUE),

  )

test_1class_ps <- test_1class %>%
  mutate(
    log_p = log10(p) %>%
      replace(., is.infinite(.), -4)
  ) %>%
  select(-data)

write_rds(test_1class_ps, "hemisphere_model_ps.rds")

## Lefty =====

test_lefties <- bsL %>%
  group_by(feature) %>%
  nest() %>%
  left_join(true_model, join_by(feature)) %>%
  mutate(

    p = map2_dbl(data, LD1,
                 ~get_p_value(tibble(LD1 = .x$LD1), .y,
                              direction = "two_sided")[[1]],
                 .progress = TRUE),

  )

test_lefties_ps <- test_lefties %>%
  mutate(
    log_p = log10(p) %>%
      replace(., is.infinite(.), -4)
  ) %>%
  select(-data)

write_rds(test_lefties_ps, "hemisphere_lefies_ps.rds")

## 4-class =====

scnd_model2 <- second_model %>%
  pivot_longer(-feature, values_to = "true_value")

test_4class <- bs2 %>%
  select(-n) %>%
  pivot_longer(-feature) %>%
  group_by(feature, name) %>%
  nest() %>%
  left_join(scnd_model2, join_by(feature, name)) %>%
  mutate(

    p = map2_dbl(data, true_value,
                  ~get_p_value(tibble(value = .x$value), .y,
                               direction = "two_sided")[[1]],
                  .progress = TRUE),

  )

qsave(test_4class, "4class_model_with_data.rds", nthreads = 8)

test_4class_ps <- test_4class %>%
  mutate(
    log_p = log10(p) %>%
      replace(., is.infinite(.), -4)
  ) %>%
  select(-data)

write_rds(test_4class_ps, "4class_model_ps.rds")
