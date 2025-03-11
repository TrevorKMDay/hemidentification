
library(MASS)
library(mltools)
library(tidyverse)
library(qs)

select <- dplyr::select

setwd("~/Google Drive/My Drive/Projects/hemisphere_fingerprinting/")

# Load data ====

all_data3 <- qread("hemiconnectome.rds")

# Create test splits ===

splits <- all_data3 %>%
  select(
    sub, handedness
  ) %>%
  distinct() %>%
  group_by(handedness) %>%
  mutate(
    group = cut_interval(row_number(), n = 5, labels = LETTERS[1:5])
  )

# Do training ====

train_on_group_notX <- function(data, test_group, outcome = "class") {

  training <- data %>%
    filter(
      group != test_group
    ) %>%
    ungroup() %>%
    select(-sub, -hemi, -handedness, -group)

  f <- as.formula(paste(outcome, "~ ."))

  lda1 <- lda(f, data = training)

  return(lda1)

}

test_on_group_X <- function(data, model, test_group) {

  test <- data %>%
    filter(
      group == test_group
    ) %>%
    ungroup() %>%
    select(-sub, -hemi, -handedness, -group)

  results <- predict(model, test)

  return(results)

}

create_contingency_table <- function(result, ground_truth) {

  tibble(truth = ground_truth, pred = result$class) %>%
    table() %>%
    as.data.frame.matrix() %>%
    return()

}

metric <- function(contingency, metric) {

  x <- as.matrix(contingency)

  if (metric == "accuracy") {

    result <- sum(diag(x)) / sum(x)

  } else if (metric == "lefty_accuracy") {

    correct <- x["leftyLH", "leftyLH"] + x["leftyRH", "leftyRH"]
    all_lefty <- sum(x[c("leftyLH", "leftyRH"), ])

    result <- correct / all_lefty

  } else if (metric == "righty_accuracy") {

    correct <- x["rightyLH", "rightyLH"] + x["rightyRH", "rightyRH"]
    all_righty <- sum(x[c("rightyLH", "rightyRH"), ])

    result <- correct / all_righty

  } else if (metric == "LH_accuracy") {

    correct <- x["leftyLH", "leftyLH"] + x["rightyLH", "rightyLH"]
    all_LH <- sum(x[c("rightyLH", "leftyLH"), ])

    result <- correct / all_LH

  } else if (metric == "RH_accuracy") {

    correct <- x["leftyRH", "leftyRH"] + x["rightyRH", "rightyRH"]
    all_RH <- sum(x[c("rightyRH", "leftyRH"), ])

    result <- correct / all_RH

  }

  return(result)

}

# TO DO: Move me into load function above
all_data4 <- splits %>%
  left_join(all_data3, by = join_by(sub, handedness)) %>%
  select(-EHI) %>%
  ungroup() %>%
  mutate(
    # Fisher Z transform correlation values - approximate normality
    across(where(is.numeric), DescTools::FisherZ)
  )

# for (i in 1:ncol(all_data4)) {
#
#   col <- unlist(all_data4[, i])
#   nans <- sum(is.nan(col))
#   if (nans > 0) message(colnames(all_data4)[i])
#
# }

write_rds(all_data4, "all_data_Z.rds")
reticulate::py_save_object(all_data4, "hemiconnectome.pickle")

for (test_group in LETTERS[1:5]) {

  # As the models get bigger, separate each one into its own call for stability
  message(paste("Starting "), test_group)
  message(date())

  model <- train_on_group_notX(all_data4, test_group)
  write_rds(model, paste0("testgroup-", test_group, "_model-lda.rds"))

  message(paste("Done training test group", test_group))
  message(date())

}

gender_model <- train_on_group_notX(all_data4, "A", outcome = "gender")
gender_pred <- test_on_group_X(all_data4, gender_model, "A")

gender_results <- create_contingency_table(gender_pred,
                                           all_data4$gender[all_data4$group == "A"])

mltools::mcc(confusionM = gender_results)

models <- tibble(test_group = LETTERS[1:5]) %>%
  mutate(
    lda = map(test_group,
                ~read_rds(paste0("testgroup-", .x, "_model-lda.rds")),
              .progress = TRUE)
  )

# Predict unseen hemispheres
results <- models %>%
  mutate(
    result = map2(lda, test_group, ~test_on_group_X(all_data4, .x, .y))
  )

contingency <- results %>%
  select(-lda) %>%
  mutate(
    table = map2(result, test_group,
                 ~create_contingency_table(.x,
                                           all_data4$class[all_data4$group == .y])),

    acc = map_dbl(table, ~metric(.x, "accuracy")),
    mcc = map_dbl(table, ~mcc(confusionM = .x)),

    lefty_acc = map_dbl(table, ~metric(.x, "lefty_accuracy")),
    righty_acc = map_dbl(table, ~metric(.x, "righty_accuracy")),
    lh_acc = map_dbl(table, ~metric(.x, "LH_accuracy")),
    rh_acc = map_dbl(table, ~metric(.x, "RH_accuracy")),


  )

metrics <- contingency %>%
  select(-result, -table) %>%
  pivot_longer(-test_group, names_to = "metric")

metric_means <- metrics %>%
  group_by(metric) %>%
  summarize(
    n = n(),
    mean = mean(value),
    sd = sd(value)
  ) %>%
  mutate(
    se = sd / sqrt(n)
  )

ggplot(metrics, aes(x = metric, y = value)) +
  geom_hline(yintercept = 1 / 4, color = "red") +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  geom_jitter(aes(color = test_group), width = 0.12, height = 0, alpha = 0.85) +
  geom_pointrange(data = metric_means,
                  aes(y = mean, ymin = mean - se, ymax = mean + se)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_discrete(labels = c("Overall accuracy", "Lefty accuracy",
                              "LH accuracy", "MCC*", "RH accuracy",
                              "Righty accuracy")) +
  theme_bw() +
  labs(x = "Metric", y = "Value", color = "Group",
       caption = "* MCC is scaled [-1, 1], where 0 represents chance assignment.")


# Differences between groups ====

t.test(contingency$lefty_acc, contingency$righty_acc, paired = TRUE)
t.test(contingency$lh_acc, contingency$rh_acc, paired = TRUE)

# Contingency Average ====

contingency_tables <- lapply(contingency$table, function(x) x / sum(x))
contingency_sum <- Reduce("+", contingency_tables)

contingency_avg <- contingency_sum / length(contingency_tables)

c_avg_plot <- contingency_avg %>%
  rownames_to_column("truth") %>%
  pivot_longer(-truth, names_to = "pred")

ggplot(c_avg_plot, aes(x = truth, y = pred, fill = value)) +
  geom_tile() +
  geom_text(aes(label = round(value, 2)), size = 10) +
  scale_fill_gradient2() +
  theme_bw() +
  labs(x = "Truth", y = "Predicted class", fill = "Rate")

# One model ====

results_A <- tibble(
    class = results$result[[1]]$class
  ) %>%
  bind_cols(results$result[[1]]$posterior) %>%
  bind_cols(results$result[[1]]$x) %>%
  bind_cols(
    all_data4 %>%
      filter(group == "A") %>%
      select(sub, handedness, hemi, class) %>%
      rename(class_true = class)
  )

ggplot(results_A, aes(x = LD1, y = LD2, color = class_true)) +
  geom_point() +
  theme_bw()

# Weights ====

Ax <- models$lda[[1]]$scaling %>%
  as_tibble(rownames = "connection") %>%
  separate_wider_delim(connection, delim = "_", names = c("ROI1", "ROI2")) %>%
  mutate(
    across(c(ROI1, ROI2), ~str_remove(., "`"))
  ) %>%
  arrange(ROI1, ROI2)

all_rois <- ordered(unique(c(Ax$ROI1, Ax$ROI2)))


Ax2 <- Ax %>%
  select(ROI1, ROI2, LD1) %>%
  mutate(

    across(starts_with("ROI"), ~ordered(.x, levels = all_rois)),
    ROIx = if_else(ROI1 > ROI2, ROI1, ROI2),
    ROIy = if_else(ROI1 > ROI2, ROI2, ROI1),

    LD1_Z = scale(LD1)[, 1],
    LD1_pepper = case_when(LD1_Z > 3 ~ 1, LD1_Z < -3 ~ -1, .default = NA)

  )

mean(Ax$LD1)
sd(Ax$LD1)
max_abs_val <- ceiling(max(abs(range(Ax2$LD1))) * 100) / 100

ggplot(Ax2, aes(LD1)) +
  geom_histogram()

ggplot(Ax2, aes(x = ROIx, y = ROIy)) +
  geom_tile(aes(fill = LD1)) +
  scale_fill_gradient2(limits = c(-max_abs_val, max_abs_val)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(Ax2, aes(x = ROIx, y = ROIy)) +
  geom_tile(aes(fill = as_factor(LD1_pepper))) +
  scale_fill_manual(values = c("red", "blue"), na.value = NA) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

library(igraph)

Ax_graph <- Ax2 %>%
  select(ROIx, ROIy, LD1) %>%
  graph_from_data_frame(directed = FALSE)

Ax_graph_3SD <- Ax2 %>%
  filter(
    abs(LD1_Z) > 3.5
  ) %>%
  select(ROIx, ROIy, LD1) %>%
  graph_from_data_frame(directed = FALSE)

plot(Ax_graph_3SD, layout = layout_nicely, edge_color = Ax_graph_3SD$LD1)

degree_centrality <- strength(Ax_graph, weights = abs(Ax2$LD1)) %>%
  as_tibble(rownames = "ROI") %>%
  arrange(desc(value))

ggplot(degree_centrality, aes(x = ROI, y = value)) +
  geom_point() +
  geom_line() +
  scale_x_discrete(limits = degree_centrality$ROI) +
  coord_flip() +
  theme_bw()

# To cifti ----

rois <- read_table("data/pscalar_info.txt", col_names = FALSE,
                   show_col_types = FALSE) %>%
  select(X3) %>%
  mutate(
    X3 = str_remove_all(X3, "^[RL]_|_ROI")
  )

out <- left_join(rois, degree_centrality, join_by(X3 == ROI))

out %>%
  select(value) %>%
  write.table("degree_centrality_to_pscalar.txt", col.names = FALSE,
              row.names = FALSE)
