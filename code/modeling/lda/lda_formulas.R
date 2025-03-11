metrics <- function(results) {

  contingency <- table(results$true_hemi, results$class)

  accuracy <- sum(diag(contingency)) / sum(contingency)
  mcc_val <- mcc(confusionM = as.matrix.data.frame(contingency))

  return(
    list("accuracy" = accuracy, "mcc" = mcc_val)
  )

}

results_to_tibble <- function(results, src) {

  results_tibble <- tibble(class = results$class) %>%
    bind_cols(results$posterior) %>%
    bind_cols(results$x) %>%
    add_column(true_hemi = src$hemi, sub = src$sub, .before = 1)

  return(results_tibble)

}
