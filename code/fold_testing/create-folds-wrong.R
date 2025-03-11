# create_split <- function(fams) {
#
#   # This code creates one fold at random
#
#   n1 <- nrow(fams)
#   ntrain <- round(0.8 * n1)
#
#   row <- 1:n1
#   training <- sample.int(n1, ntrain)
#   test <- row[!(row %in% training)]
#
#   train_set <- unnest(fams[training, ], data)
#   test_set <- unnest(fams[test, ], data)
#
#   return(tibble(train = list(train_set), test = list(test_set)))
#
# }



create_folds <- function(fams) {

  # We proposed to create five *mutually exclusive* folds. This function will
  #   do that. It repeats the process if the resulting sizes are unevenly
  #   distributed according to a chi^2 test.

  conds <- LETTERS[1:5]
  equal_split <- nrow(fams) / length(conds)

  # loop control vars
  equally_split <- FALSE

  while (!equally_split) {

    fams$fold <- cluster_ra(fams$family, conditions = LETTERS[1:5])

    fams2 <- fams %>%
      group_by(fold) %>%
      nest() %>%
      mutate(
        size = map_int(data, nrow),
        mean_age = map_dbl(data, ~mean(.x$age)),
        hand1R_ratio = map_dbl(data, ~sum(.x$hand1 == "right")) / size,
        hand2R_ratio = map_dbl(data, ~sum(.x$hand2 == "right2")) / size,
      ) %>%
      arrange(fold)

    # Check that the sizes are equally distributed, because the families are
    #   different sizes
    size_test <- chisq.test(tibble(sizes = fams2$size, equal = equal_split),
                            simulate.p.value = TRUE)

    age_test <- chisq.test(tibble(ages = fams2$mean_age,
                                  equal = mean(fams$age)),
                           simulate.p.value = TRUE)

    hand1R_test <- chisq.test(
      tibble(hand1R = fams2$hand1R_ratio,
             equal = sum(fams$hand1 == "right") / nrow(fams)),
      simulate.p.value = TRUE
    )

    hand2R_test <- chisq.test(
      tibble(hand2R = fams2$hand2R_ratio,
             equal = sum(fams$hand2 == "right2") / nrow(fams)),
      simulate.p.value = TRUE
    )

    # Try again if the sample is unbalanced
    if (all(c(size_test$p.value, age_test$p.value, hand1R_test$p.value,
              hand2R_test$p.value) > 0.05))
      equally_split <- TRUE

  }

  fams3 <- select(fams2, fold, data)

  return(fams3)

}

split_testing <- create_folds(rdat2)

create_tt_splits <- function(fold) {

  fold$split <- block_ra(blocks = fold$hand1, conditions = c("train", "test"),
                         prob_each = c(0.8, 0.2))

  split_sizes <- nrow(fold) * c(0.8, 0.2)

  # loop control vars
  equally_split <- FALSE

  while (!equally_split) {

    fold2 <- fold %>%
      group_by(split) %>%
      nest()%>%
      mutate(
        size = map_int(data, nrow),
        mean_age = map_dbl(data, ~mean(.x$age)),
        hand1R_ratio = map_dbl(data, ~sum(.x$hand1 == "right")) / size,
        hand2R_ratio = map_dbl(data, ~sum(.x$hand2 == "right2")) / size,
      ) %>%
      arrange(split)

    # Check that the sizes are equally distributed, because the families are
    #   different sizes
    size_test <- chisq.test(tibble(sizes = fold2$size, equal = split_sizes),
                            simulate.p.value = TRUE)

    age_test <- chisq.test(tibble(ages = fold2$mean_age,
                                  equal = mean(fold$age)),
                           simulate.p.value = TRUE)

    hand1R_test <- chisq.test(
      tibble(hand1R = fold2$hand1R_ratio,
             equal = sum(fold$hand1 == "right") / nrow(fold)),
      simulate.p.value = TRUE
    )

    hand2R_test <- chisq.test(
      tibble(hand2R = fold2$hand2R_ratio,
             equal = sum(fold$hand2 == "right2") / nrow(fold)),
      simulate.p.value = TRUE
    )

    # Try again if the sample is unbalanced
    if (all(c(size_test$p.value, age_test$p.value, hand1R_test$p.value,
              hand2R_test$p.value) > 0.05))
      equally_split <- TRUE

  }

  fold3 <- select(fold2, split, data)

  return(fold3)

}

splittt_testing <- split_testing %>%
  mutate(
    train_test = map(data, ~create_tt_splits(.x))
  ) %>%
  select(-data) %>%
  unnest(train_test) %>%
  unnest(data)

