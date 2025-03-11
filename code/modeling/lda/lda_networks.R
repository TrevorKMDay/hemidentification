library(tidyverse)
library(MASS)
library(mltools)
library(ggrepel)

select <- dplyr::select

setwd("/Volumes/thufir/HCPYA")

# Setup ====

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
                    sep = "_"),
    connection_type = if_else(connection_type == "ips_withinnet",
                              paste0(connection_type, "_", net1),
                              connection_type)
  ) %>%
  select(roi_roi, net1, net2, connection_type) %>%
  distinct() %>%
  ungroup() %>%
  na.omit()

# List of all networks
networks <- unique(ctypes$net1)

training_set <- read_rds("hemiconns/A_train_RH_n30_ALL8_wide.rds")
test_set <- read_rds("hemiconns/A_test_RH_n30_ALL8_wide.rds")

source("code/lda_formulas.R")

# Embed training sets ====

column_selector <- tibble(cnctn = c(paste0("ips_withinnet_", networks),
                                    "ips_betweennet"))

training_sets <-  column_selector %>%
  mutate(
    training_data = map(cnctn,
                        ~select(training_set,
                                any_of(c("sub", "hemi",
                                        ctypes$roi_roi[ctypes$connection_type == .x]))
                                )
                        ),
    n_cols = map_int(training_data, ncol)
  )

training_models <- training_sets %>%
  filter(
    # This keeps all networks larger than and including Aud.
    n_cols >= 30
  )  %>%
  mutate(
    lda = map(training_data, ~lda(hemi ~ ., data = select(.x, -sub)),
              .progress = TRUE)
  )

test_sets <- column_selector %>%
  mutate(
    test_data = map(cnctn,
                      ~select(test_set,
                              any_of(c("sub", "hemi",
                                       ctypes$roi_roi[ctypes$connection_type == .x]))
                      )
                    ),
    n_cols = map_int(test_data, ncol)
  ) %>%
  filter(
    # This keeps all networks larger than and including Aud.
    n_cols >= 30
  ) %>%
  left_join(
    select(training_models, -training_data),
    join_by(cnctn, n_cols)
  ) %>%
  mutate(
    pred = map2(lda, test_data, ~predict(.x, .y)),
    results = map2(pred, test_data, ~results_to_tibble(.x, .y)),
    metrics = map(results, metrics),
    accuracy = map_dbl(metrics, ~ .x$accuracy),
    mcc = map_dbl(metrics, ~ .x$mcc)
  )

test_sets_metrics <- test_sets %>%
  select(cnctn, n_cols, accuracy, mcc) %>%
  pivot_longer(c(accuracy, mcc), names_to = "metric") %>%
  mutate(
    cnctn = str_remove_all(cnctn, "ips_withinnet_|ips_")
  )

ggplot(test_sets_metrics, aes(x = n_cols, y = value)) +
  geom_point(aes(color = cnctn)) +
  geom_smooth(group = 1, method = "lm") +
  geom_label_repel(aes(label = cnctn, fill = cnctn)) +
  scale_x_log10() +
  coord_cartesian(ylim = c(0, 1)) +
  scale_fill_discrete(guide = "none") +
  facet_wrap(vars(metric)) +
  theme_bw() +
  labs(x = "# variables", y = "Metric value", color = "Network")

bootstrap <- function(train, test, size) {

  cols <- colnames(train)[-c(1, 2)]
  random_cols <- sample(cols, size = size)

  temp_training <- train %>%
    select(hemi, all_of(random_cols))

  temp_test <- test  %>%
    select(hemi, all_of(random_cols))

  model <- lda(hemi ~ ., temp_training)
  results <- results_to_tibble(predict(model, temp_test),
                               src = temp_training)

  metrics_vals <- as_tibble(metrics(results))

  return(metrics_vals)

}

lang_bootstraps <- tibble(rep = 1:20) %>%
  mutate(
    results = map(rep, ~bootstrap(train = training_set, test = test_set,
                                  size = 93))
  ) %>%
  unnest(results) %>%
  pivot_longer(-rep, names_to = "metric")

random_assignment <- tribble(
  ~metric, ~value,
  "accuracy", 0.5,
  "mcc", 0
)

p1 <- ggplot(NULL) +
  geom_hline(data = random_assignment, aes(yintercept = value),
             color = "red") +
  geom_point(data = test_sets_metrics,
             aes(x = n_cols, y = value, color = cnctn)) +
  geom_label_repel(data = test_sets_metrics,
                   aes(x = n_cols, y = value, label = cnctn, fill = cnctn)) +
  scale_x_log10(limits = c(10, NA)) +
  coord_cartesian(ylim = c(NA, 1)) +
  scale_fill_discrete(guide = "none") +
  scale_color_discrete(guide = "none") +
  facet_wrap(vars(metric), scales = "free_y") +
  theme_bw() +
  labs(x = "# variables", y = "Metric value", color = "Network")

p1 +
  geom_jitter(
    data = lang_bootstraps, aes(x = 93, y = value),
    width = 0.05, height = 0, alpha = 0.25
  )
