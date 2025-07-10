library(tidyverse)

memory <- read_tsv("../data/mem_usage.tsv", col_names = c("timestamp", "usage"),
                   show_col_types = FALSE) %>%
  mutate(
    timestamp = timestamp %>%
      str_remove(" EDT 2025") %>%
      as_datetime(format = "%a %b %d %I:%M:%S %p"),
    usage_kb = as.numeric(str_remove(usage, "kB")),
    usage_mb = usage_kb / 1000
  )

ggplot(memory, aes(x = timestamp, y = usage_mb)) +
  geom_line()

memory2 <- memory %>%
  filter(
    timestamp < as_datetime("2025-03-20 10:30:00")
  ) %>%
  mutate(
  )


