library(tidyverse)

setwd("/Users/tkmd/MyDrive/Projects/hemisphere_fingerprinting/code/modeling_new/05-flipped_parcels/")
source("../../bootstrap_mcc.R")

library(parallel)
library(doParallel)
library(foreach)

cores <- detectCores()
cl <- makeCluster(cores[1] - 1) #not to overload your computer
registerDoParallel(cl)

files <- list.files("results/", full.names = TRUE)

data0 <- tibble(f = files) %>%
  mutate(
    data = map(f, read_csv, show_col_types = FALSE)
  )

data <- data0 %>%
  mutate(
    flip = str_extract(f, "fp.."),
    hand = str_extract(f, "[LR]HD")
  )  %>% 
  select(flip, hand, data) %>%
  mutate(

    mcc = map(data, bootstrap_mcc, ground = "hemi")

  )

results_mcc <- data %>%
  select(-data) %>%
  unnest(mcc)

ggplot(results_mcc, aes(x = flip, y = mcc_mean)) +
  geom_pointrange(
    aes(ymin = mcc_ci_lo, ymax = mcc_ci_hi, color = hand),
    position = position_dodge2(width = 0.1)
  ) +
  theme_bw()
