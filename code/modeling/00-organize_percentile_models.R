library(tidyverse)
library(qs)
library(reticulate)

select <- dplyr::select

setwd("~/Google Drive/My Drive/Projects/hemisphere_fingerprinting/code/modeling/")

hc <- qread("inputs/hemiconnectome.rds")

connections <- read_rds("../../data/residualized_hemiconnectome.rds")

connections1 <- connections %>%
  mutate(
    conn = paste(ROI1, ROI2, sep = "_")
  ) %>%
  select(conn, rZ_rzd) %>%
  mutate(
    pctile = percent_rank(rZ_rzd)
  )

pctiles <- c(seq(50, 90, by = 5), 95, 97.5, 99, 99.5, 99.75, 99.9) / 100

dir.create("inputs/subsets/", showWarnings = FALSE, recursive = TRUE)

rev_conn <- function(x) {

  y <- str_split_1(x, "_")

  return(paste(y[2], y[1], sep = "_"))

}

for (p in pctiles) {

  temp <- connections1 %>%
    filter(
      pctile >= p
    )

  conns <- temp$conn
  rev_conns <- Vectorize(rev_conn)(conns)

  new_hc <- hc %>%
    select(sub, group, EHI, gender, age_group, hemi, starts_with("handedness"),
           starts_with("class"), any_of(conns), any_of(rev_conns))

  name <- paste0("hemiconnectome_pctile", p * 100, ".pickle")

  py_save_object(new_hc, paste0("inputs/subsets/", name))

  message(name)

}
