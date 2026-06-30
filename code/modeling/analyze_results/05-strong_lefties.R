library(tidyverse)
library(mltools)

setwd("~/Google Drive/My Drive/Projects/hemisphere_fingerprinting/code/modeling/analyze_results/")

strong_lefty_files <- list.files("../results/strong_lefties/",
                                 pattern = "_test-[0-9]*_.*.csv",
                                 full.names = TRUE)

results0 <- tibble(f = strong_lefty_files) %>%
  mutate(
    data = map(f, read_csv, show_col_types = FALSE)
  ) %>%
  unnest(data)

results <- results0 %>%
  select(-`...1`)

conf <- table(results$ground, results$predicted)

mcc(confusionM = as.matrix.data.frame(conf))
