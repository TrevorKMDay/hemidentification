
ehi_files <- list.files("results/ehi/", pattern = "*.csv",
                        full.names = TRUE)

results <- tibble(f = ehi_files) %>%
  mutate(

    data = map(f, read_csv, show_col_types = FALSE),
    f = basename(f),

    method = str_extract(f, "method-[^_]*") %>%
      str_remove("method-"),

    # outcome = str_extract(f, "outcome-[^_]*"),

    group = str_extract(f, "test-[^_]*") %>%
      str_remove("test-"),

    input = str_extract(f, "data-.*") %>%
      str_remove_all("data-|.csv"),

  ) %>%
  select(-f) %>%
  unnest(data) %>%
  select(-`...1`)

ggplot(results, aes(x = ground, y = predicted)) +
  geom_jitter(aes(color = group),
              height = 0, width = 1,
              alpha = 0.5) +
  geom_smooth(method = "lm", color = "black") +
  geom_hline(yintercept = 0, color = "red",
             linetype = "solid") +
  geom_hline(yintercept = 100, color = "red",
             linetype = "dashed") +
  scale_x_continuous(breaks = seq(-100, 100, 50)) +
  scale_y_continuous(limits = c(-100, NA)) +
  theme_bw() +
  facet_grid(rows = vars(input),
             cols = vars(method))
