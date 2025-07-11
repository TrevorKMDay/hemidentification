library(tidyverse)
library(qs)

# This seems to take ~1 h per 5,000 reps (so this whole thing could take 4 h)
# with 2x10,000 reps, so separating it out

bs_file <- "lda_full_bootstrapped.rds"

if (!file.exists(bs_file)) {

  bs_f <- list.files("../modeling/lda_full/", "full_bs[0-9]*_scalings.csv",
                     full.names = TRUE)

  bs0 <- tibble(f = bs_f) %>%
    mutate(
      data = map(f, read_csv, show_col_types = FALSE, .progress = TRUE),
      n = row_number()
    )

  bs <- bs0 %>%
    select(n, data) %>%
    unnest(data)

  qsave(bs, bs_file)

} else {

  bs <- qread(bs_file, nthreads = 2)

}

### 4-class model

bs_file2 <- "lda_full_bootstrapped2.rds"

if (!file.exists(bs_file2)) {

  bs_fs2 <- list.files("../modeling/lda_full/",
                       "full_class_bs[0-9]*_scalings.csv",
                       full.names = TRUE)

  bs2_0 <- tibble(f = bs_fs2) %>%
    mutate(
      data = map(f, read_csv, show_col_types = FALSE,
                 .progress = TRUE),
      n = row_number()
    )

  bs2 <- bs2_0 %>%
    select(n, data) %>%
    unnest(data)

  qsave(bs2, bs_file2)

} else {

  bs2 <- qread(bs_file2, nthreads = 2)

}
