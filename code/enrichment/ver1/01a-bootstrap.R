library(tidyverse)
library(qs)
library(data.table)

setwd("~/MyDrive/Projects/hemisphere_fingerprinting/code/enrichment/")

# This seems to take ~1 h per 5,000 reps (so this whole thing could take 4 h)
# with 2x10,000 reps, so separating it out

bs_file <- "lda_full_bootstrapped.rds"

# one_class_glob <- "../modeling/lda_full/full_bs0[0-9]*_scalings.csv"
#
# x <- fread(cmd = paste("awk 'NR==1||FNR!=1'", one_class_glob),
#             skip = 0, col.names = c("feature", "LD1"), header = TRUE,
#             verbose = TRUE, showProgress = TRUE)
#
# which(is.na(as.numeric(x$LD1)))

if (!file.exists(bs_file)) {

  bs_f <- list.files("../modeling/lda_full/bootstrapping/",
                     "full_bs[0-9]*_scalings.csv",
                     full.names = TRUE)

  # bs0 <- tibble(f = bs_f) %>%
  #   mutate(
  #     data = map(f, read_csv, show_col_types = FALSE, .progress = TRUE),
  #     n = row_number()
  #   )

  bs0 <- tibble(f = bs_f) %>%
    mutate(
      data = map(f, fread, .progress = TRUE),
      n = row_number()
    )

  bs <- bs0 %>%
    select(n, data) %>%
    unnest(data)

  qsave(bs, bs_file)

} else {

  bs <- qread(bs_file, nthreads = 2)

}


### Lefty model

bs_fileL <- "lda_full_lefies.rds"

if (!file.exists(bs_fileL)) {

  bs_fs2 <- list.files("../modeling/lda_full/",
                       "full_lefties_bs[0-9]*_scalings.csv",
                       full.names = TRUE)

  bs2_0 <- tibble(f = bs_fs2) %>%
    mutate(
      data = map(f, fread, .progress = TRUE),
      n = row_number()
    )

  bs2 <- bs2_0 %>%
    select(n, data) %>%
    unnest(data)

  qsave(bs2, bs_fileL)

} else {

  bs2 <- qread(bs_file2, nthreads = 2)

}

ids <- str_extract(bs_fs2, "bs[0-9]*") %>%
  str_remove("bs")

true_ids <- str_pad(0:9999, 4, "left", "0")

setdiff(true_ids, ids)


### 4-class model

bs_file2 <- "lda_full_bootstrapped2.rds"

if (!file.exists(bs_file2)) {

  bs_fs2 <- list.files("../modeling/lda_full/",
                       "full_class_bs[0-9]*_scalings.csv",
                       full.names = TRUE)

  bs2_0 <- tibble(f = bs_fs2) %>%
    mutate(
      data = map(f, fread, .progress = TRUE),
      n = row_number()
    )

  bs2 <- bs2_0 %>%
    select(n, data) %>%
    unnest(data)

  qsave(bs2, bs_file2)

} else {

  bs2 <- qread(bs_file2, nthreads = 2)

}

ids <- str_extract(bs_fs2, "bs[0-9]*") %>%
  str_remove("bs")

true_ids <- str_pad(0:9999, 4, "left", "0")

setdiff(true_ids, ids)
