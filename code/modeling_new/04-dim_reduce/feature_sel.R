library(tidyverse)
library(reticulate)

# Setup and get data =====

setwd(str_glue("/Users/tkmd/MyDrive/Projects/hemisphere_fingerprinting/code/",
               "modeling_new/04-dim_reduce"))

results <- py_load_object("feature_sel_results.pickle")

format_result <- function(x) {

  # Extract the data from the Python object 

  name <- x[[1]][[1]]
  size <- x[[2]]

  group_mccs <- x[[3]][[2]]
  overall_mcc <- x[[3]][[3]]

  mcc <- bind_rows(group_mccs, overall_mcc) %>%
    as_tibble() %>%
    mutate(
        name = name,
        k = size,
        group = replace_na(group, "overall")
    )
  
}

results_fmt <- map(results, format_result) %>%
  bind_rows() %>%
  mutate(
    k_base = log(k, base = 2)
  )

# Break out fold and overall results
results_fmt_groups <- filter(results_fmt, group != "overall")
results_fmt_overall <- filter(results_fmt, group == "overall")

rule_of_thumbs <- tribble(
    ~name, ~rule10x, ~rulesqrt,
    "hclhd", 18, 14,
    "hcrhd", 168, 42
  ) %>%
  pivot_longer(-name, names_to = "rule") %>%
  mutate(
    value = log(value, base = 2)
  )

ggplot(results_fmt_groups, aes(x = k_base, y = mcc, color = name)) +
  geom_vline(
    data = rule_of_thumbs, 
    aes(xintercept = value, color = name, linetype = rule),
    linewidth = 1
  ) +
  geom_hline(
    yintercept = 0, 
    color = "red"
  ) +
  geom_point(
    # Background points
    alpha = .5,
    shape = 18
  ) +
  geom_pointrange(
    data = results_fmt_overall, 
    aes(ymin = mcc_ci_lo, ymax = mcc_ci_hi),
  ) +
  geom_line(
    data = results_fmt_overall, 
  ) +
  scale_x_continuous(
    breaks = seq(0, 12, by = 3),
    labels = 2 ** seq(0, 12, by = 3),
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2)
  ) +
  scale_color_discrete(
    labels = c("hclhd" = "Sinistral", "hcrhd" = "Dextral")
  ) +
  scale_linetype(
    labels = c("rule10x" = "10x", "rulesqrt" = "sqrt")
  ) +
  theme_bw() +
  labs(x = "# features (base-2 log-transformed)", y = "MCC (95% CI)",
        color = "Handedness", linetype = "# features") +
  theme(legend.position = "bottom")

ggsave("plots/feature_selection_LR.png", width = 6.5, height = 4)

# Features

lefty_features <- read_csv("lefty_scores.csv") 
righty_features <- read_csv("righty_scores.csv")

features_all <- left_join(lefty_features, righty_features) %>%
  mutate(
    # Get feature name and order from largest to smallest
    feature = str_remove(feature, "^cx_"),
    lefty_rank = rank(-lefty_score),
    righty_rank = rank(-righty_score)
  ) 

# Top 16 papers for manuscript

features <- features_all %>%
  filter(
    lefty_rank <= 16 | righty_rank <= 16
  )

features2 <- features %>%
  select(-ends_with("score")) %>%
  arrange(righty_rank) %>%
  separate_wider_delim(feature, delim = "_", names = c("roi1", "roi2"),
                        cols_remove = FALSE)

clipr::write_clip(features2)

# Which regions appear more than once?
table(c(features2$roi1, features2$roi2))

# Get the symmetry scores and rename it to symmetry to be clearer
overlap_index <- read_csv("~/MyDrive/Projects/parcel_overlap/glasser_LR_overlap.csv",
                            show_col_types = FALSE) %>%
  mutate(
    avg_size = (sizeL + sizeR) / 2
  ) %>%
  select(label, dice, avg_size) %>%
  rename(
    symmetry = dice
  )

# Get each endpoint's asymmetry score and put it together 
features3 <- features_all %>%
  separate_wider_delim(feature, delim = "_", names = c("roi1", "roi2")) %>%
  left_join(
    overlap_index, by = join_by(roi1 == label)
  ) %>%
  left_join(
    overlap_index, by = join_by(roi2 == label), suffix = c("1", "2")
  ) %>%
  mutate(

    total_symmetry = (symmetry1 + symmetry2) / 2,
    total_size = (avg_size1 + avg_size2) / 2,
    
    # Z transform because the sinistral model is on different scale than 
    # righties
    across(c(ends_with("score"), starts_with("total")), ~scale(.x)[,1], 
            .names = "{col}Z")

  )

lm_r <- lm(righty_scoreZ ~ total_symmetryZ + total_sizeZ, data = features3)
lm_l <- lm(lefty_scoreZ ~ total_symmetryZ + total_sizeZ, data = features3)

features4 <- features3 %>%
  select(roi1, roi2, ends_with("scoreZ"), starts_with("total")) %>%
  pivot_longer(ends_with("scoreZ")) %>%
  mutate(
    name = replace_values(
      name,
      "lefty_scoreZ" ~ "Sinistrals",
      "righty_scoreZ" ~ "Dextrals"
    )
  )

features5 <- features3 %>%
  select(roi1, roi2, ends_with("rank"), starts_with("total")) %>%
  pivot_longer(ends_with("rank"))

ggplot(features4, aes(x = total_symmetryZ, y = value)) +
  geom_point(alpha = 0.2, aes(color = total_sizeZ)) +
  geom_smooth(method = "lm") +
  viridis::scale_color_viridis() +
  facet_wrap(vars(name), scales = "free_y") +
  theme_bw() +
  labs(x = "Z(Combined symmetry)", y = "Z(F classifier score)", color = "Z(Total size)")

ggsave("plots/fclassifierscore_asymmetry.png", width = 6.5, height = 3)

ggplot(features5, aes(x = total_symmetryZ, y = value)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm") +
  facet_wrap(vars(name), scales = "free_y") +
  theme_bw()

features3_toprank <- features3 %>%
  filter(
    righty_rank <= 16 | lefty_rank <= 16
  )
