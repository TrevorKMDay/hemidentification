# Breaking out some follow-ups to 00-organize_data.R for parsimony

library(tidyverse)
library(qs)
library(reticulate)

select <- dplyr::select

setwd("~/Google Drive/My Drive/Projects/hemisphere_fingerprinting/code/modeling/")

hc1 <- qread("inputs/hemiconnectome.rds")

# Check these data - if range > [-1, 1], then it's FisherZ transformed
# hc1 %>%
#   select(where(is.numeric), -EHI) %>%
#   range()

# Connection normalization ====

# Normalize each value by connection, i.e., normalize sub-1 connection #1
# against all other connection #1s.
# Normalizes Fisher-Z values

# Reorder columns and remove "2" cols, not really using that naymore
hc2 <- hc1 %>%
  select(sub, group, gender, age_group, EHI, handedness, class, hemi,
          contains("_")) %>%
  mutate(
    # Normalize between columns
    across(c(where(is.numeric), -EHI), ~scale(.x)[,1])
  )


qsave(hc2, "inputs/hemiconnectome_connNormal.rds")
py_save_object(hc2, "inputs/hemiconnectome_connNormal.pickle")
