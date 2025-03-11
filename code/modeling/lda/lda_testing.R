library(tidyverse)
library(MASS)
library(mltools)

select <- dplyr::select

# Read in connection types
#   This contains both A_B and B_A, oh well
ctypes <- read_rds("data/glasserhcp/allpairs_connectiontype.rds") %>%
  filter(
    # Only working with ipsilateral connections, remove them here to speed
    #   things up
    laterality == "ips"
  ) %>%
  mutate(
    roi_roi = paste(str_remove(roi1, "^._"), str_remove(roi2, "^._"),
                    sep = "_")
  ) %>%
  select(roi_roi, connection_type) %>%
  distinct()

training_set <- read_csv("/Volumes/thufir/HCPYA/hemiconns/A_train_RH_n30_ALL8_wide.csv",
                          show_col_types = FALSE)

write_rds(training_set, "hemiconns/A_train_RH_n30_ALL8_wide.rds")

test_set <- read_csv("/Volumes/thufir/HCPYA/hemiconns/A_test_RH_n30_ALL8_wide.csv",
                 show_col_types = FALSE, progress = TRUE)

write_rds(test_set, "hemiconns/A_test_RH_n30_ALL8_wide.rds")



# Train the models ====

lda_input <- select(training_set, -sub)

lda_input_1each <- training_set %>%
  filter(
    row_number() <= 30 & hemi == "R" |
    row_number() >  30 & hemi == "L",
  ) %>%
  select(-sub)

lda_input_half <- training_set %>%
  group_by(hemi) %>%
  slice_sample(prop = 0.5) %>%
  select(-sub)

# Full model
hemi_lda <- lda(hemi ~ ., lda_input, tol = 1e-6)

# one hemisphere per particpant
hemi_lda_1each <- lda(hemi ~ ., lda_input_1each)

hemi_lda_half <- lda(hemi ~ ., lda_input_half)

# Evaluate results ====

self_output_all <- predict(hemi_lda, lda_input)
self_output_1e  <- predict(hemi_lda_1each, lda_input)
self_output_half <- predict(hemi_lda_half, lda_input)

results_self_all <- results_to_tibble(self_output_all, lda_input)
results_self_1e <- results_to_tibble(self_output_1e, lda_input)
results_self_half <- results_to_tibble(self_output_half, lda_input)

metrics_self_all <- metrics(results_self_all)
metrics_self_1e <- metrics(results_self_1e)
metrics_self_half <- metrics(results_self_half)

# Out-of-sample ====

test_output_all <- predict(hemi_lda, test_set)
test_output_1e  <- predict(hemi_lda_1each, test_set)
test_output_half <- predict(hemi_lda_half, test_set)

results_all <- results_to_tibble(test_output_all, test_set)
results_1e <- results_to_tibble(test_output_1e, test)
results_half <- results_to_tibble(test_output_half, test)

metrics_all <- metrics(results_all)
metrics_1e <- metrics(results_1e)
metrics_half <- metrics(results_half)

# Compare me

all_results <- tibble(results = list(metrics_all, metrics_1e, metrics_half,
                                     metrics_self_all, metrics_self_1e,
                                     metrics_self_half)) %>%
  mutate(
    size = rep(c(60, 30, 30), 2),
    sample = rep(c("test", "self"), each = 3),
    model = rep(c("all", "1e", "half"), 2),
    accuracy = map_dbl(results, ~ .x$accuracy),
    mcc = map_dbl(results, ~ .x$mcc)
  ) %>%
  select(-results) %>%
  pivot_longer(c(accuracy, mcc))

ggplot(all_results, aes(x = model, y = value, color = sample)) +
  geom_point() +
  geom_line(aes(group = interaction(sample, name),
                linetype = name)) +
  scale_x_discrete(limits = c("half", "1e", "all")) +
  scale_y_continuous(limits = c(NA, 1)) +
  theme_bw() +
  facet_wrap(vars(name), scales = "free")

