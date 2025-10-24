bootstrap_mcc <- function(data, R = 1000, ground = "ground",
                          predicted = "predicted", levels = NULL,
                          method = 1) {

  estimates <- rep(NA, R)

  estimates <- foreach (i = 1:R, .combine = "c") %dopar% {

    new_data <- dplyr::slice_sample(data, prop = 1, replace = TRUE)

    g <- unlist(new_data[, ground], use.names = FALSE)
    p <- unlist(new_data[, predicted], use.names = FALSE)

    values <- unique(c(g, p))

    if (!is.null(levels))
      if (length(values) != length(levels))
        warning("Levels do not match values in data")
    else
      levels <- values

    g <- factor(g, levels = levels)
    p <- factor(p, levels = levels)

    conf <- as.data.frame.matrix(table(g, p))
    conf2 <- conf[levels, levels]

    mltools::mcc(confusionM = conf2)

  }

  est_mean <- mean(estimates)

  if (method == 1) {

    # Preserving this method for later
    est_sd <- sd(estimates)
    est_ci <- est_sd * 2

    result <- tibble(mcc_mean = est_mean,
                     mcc_ci_lo = est_mean - est_ci,
                     mcc_ci_hi = est_mean + est_ci)

  } else if (method == 2) {

    est_ci <- unname(quantile(estimates, c(0.025, 0.975)))

    result <- tibble(mcc_mean = est_mean,
                     mcc_ci_lo = est_ci[1],
                     mcc_ci_hi = est_ci[2])

  }

  return(result)

}
