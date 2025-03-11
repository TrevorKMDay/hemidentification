library(tidyverse)
library(MASS)
library(mltools)

select <- dplyr::select

training_set <- read_csv("/Volumes/thufir/HCPYA/hemiconns/A_train_RH_n30_RS1_wide.csv",
                         show_col_types = FALSE)

test_set <- read_csv("/Volumes/thufir/HCPYA/hemiconns/A_test_RH_n30_RS1_wide.csv",
                         show_col_types = FALSE)

confusion <- function(model_results_df) {
  
  cmatrix <- table(model_results_df[, c("true_hemi", "class")]) %>%
    as_tibble() %>%
    pivot_wider(names_from = true_hemi, values_from = n) %>%
    column_to_rownames("class")
  
  return(cmatrix)
  
}

# Train the models ====

# Full model
hemi_lda <- lda(hemi ~ ., select(training_set, -sub))

model_results <- predict(hemi_lda, test_set)
model_results_df <- results_to_tibble(model_results, test_set)
model_metrics <- metrics(model_results_df)

confusion(model_results_df)
