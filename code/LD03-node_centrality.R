library(tidyverse)
library(igraph)
library(corrr)

# Setup ====

setwd("~/MyDrive/Projects/hemisphere_fingerprinting/code/")
select <- dplyr::select

# Load data ====

scaling_files <- list.files("modeling/results/base", "*_scalings.csv",
                            full.names = TRUE)

scalings <- tibble(f = scaling_files) %>%
  mutate(
    data = map(f, read_csv, show_col_types = FALSE)
  ) %>%
  unnest(data) %>%
  mutate(
    group = str_extract(f, "test-.") %>%
      str_remove("test-")
  ) %>%
  group_by(group) %>%
  mutate(
    across(starts_with("LD"), percent_rank, .names = "{.col}_pctile")
  ) %>%
  select(f)

scalings <- read_csv("modeling/lda_full/full_scalings.csv",
                     show_col_types = FALSE)

scalings_avg_pctile <- scalings %>%
  select(feature, group, ends_with("pctile")) %>%
  pivot_longer(-c(feature, group)) %>%
  group_by(feature, name) %>%
  summarize(
    n = n(),
    mean_value = mean(value),
    sd_value = sd(value)
  ) %>%
  mutate(
    se_value = sd_value / sqrt(n)
  ) %>%
  separate_wider_delim(feature, delim = "_", names = c("ROI1", "ROI2"),
                       cols_remove = FALSE) %>%
  ungroup()

scalings_avg_pctile_wide <- scalings_avg_pctile %>%
  pivot_wider(names_from = name, values_from = mean_value, id_cols = feature)

# Correlations are low, as they should be, these should be somewhat orthogonal
correlate(scalings_avg_pctile_wide)

## Load CA RSNs ====

networks0 <- read_csv(file = "../data/ColeAnticevicNetPartition/all_parcels_to_networks.csv",
                      show_col_types = FALSE)

nets <- networks0 %>%
  mutate(
    # Use LH network assignment, and combine Visual networks 1 and 2
    network = if_else(L == "Visual2", "Visual", L),
    net_short = case_match(
      network,
      "Auditory" ~ "AudN",
      "CinguloOpercular" ~ "CON",
      "Default" ~ "DMN",
      "DorsalAttention" ~ "DAN",
      "Frontoparietal" ~ "FPN",
      "Language" ~ "LN",
      "OrbitoAffective" ~ "OAN",
      "PosteriorMultimodal" ~ "pMN",
      "Somatomotor" ~ "SMN",
      "VentralMultimodal" ~ "vMN",
      "Visual" ~ "VisN"
    )
  ) %>%
  rename(
    roi = label2
  ) %>%
  select(roi, network, net_short)

nets_size <- nets %>%
  group_by(network) %>%
  summarize(
    size = n()
  )

exclude_nets <- c("OrbitoAffective", "PosteriorMultimodal", "VentralMultimodal")

## Combine. ====

scalings_net <- left_join(scalings_avg_pctile, nets,
                          join_by("ROI1" == "roi")) %>%
  left_join(nets, join_by("ROI2" == "roi"), suffix = c("1", "2")) %>%
  mutate(
    within_net = (network1 == network2)
  )

ggplot(scalings_net, aes(x = name, y = mean_value, fill = within_net)) +
  geom_boxplot() +
  theme_bw()

lm(mean_value ~ name*within_net, data = scalings_net) %>%
  summary()

## Load data ====

extract_node_strength <- function(connections) {

  graph <- graph_from_data_frame(connections, directed = FALSE) %>%
    simplify()

  nodes <- names(V(graph))
  E(graph)$weight <- abs(connections$scaling)

  strength_df <- tibble(node = nodes, strength = strength(graph)) %>%
    arrange(node) %>%
    mutate(
      z = scale(strength)[, 1],
      q = percent_rank(strength)
    )

  return(strength_df)

}

# Create matrices =====

LD1_all <- scalings %>%
  group_by(feature) %>%
  summarize(
    LD1 = mean(LD1)
  ) %>%
  separate_wider_delim(feature, delim = "_", names = c("ROI1", "ROI2")) %>%
  arrange(ROI1, ROI2) %>%
  left_join(nets, join_by("ROI1" == "roi")) %>%
  left_join(nets, join_by("ROI2" == "roi"), suffix = c("1", "2")) %>%
  mutate(
    within_net = (network1 == network2)
  )

LD1_netsonly <- LD1_all %>%
  filter(
    within_net
  )

strength1 <- LD1_all %>%
  rename(
    scaling = LD1
  ) %>%
  extract_node_strength()

strength2 <- LD1_netsonly %>%
  rename(
    scaling = LD1
  ) %>%
  extract_node_strength()

strengths <- full_join(strength1, strength2, suffix = c("1", "2"),
                       join_by(node)) %>%
  full_join(nets, join_by(node == roi)) %>%
  filter(
    !(network %in% exclude_nets)
  ) %>%
  rename(
    net = net_short
  ) %>%
  fastDummies::dummy_cols("net")

ggplot(strengths, aes(x = strength1, y = strength2, color = net)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw() +
  labs(title = "Node strength (all, nets)")

ggplot(strengths, aes(x = q1, y = q2, color = net)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw() +
  labs(title = "Node strength quantile (all, nets)")

# Models ====

lm1_all <- lm(z1 ~ net_AudN + net_CON + net_DAN + net_DMN + net_FPN + net_LN +
                net_SMN, data = strengths)

lm1_wnet <- lm(z2 ~ net_AudN + net_CON + net_DAN + net_DMN + net_FPN + net_LN +
                net_SMN, data = strengths)

strengths2 <- strengths %>%
  pivot_longer(c(strength1, strength2, z1, z2, q1, q2)) %>%
  mutate(
    model = if_else(str_extract(name, ".$") == 1, "all", "wnet"),
    name = str_remove(name, "[12]$")
  ) %>%
  pivot_wider()

ggplot(strengths2, aes(x = net, y = q, fill = model)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(NA, NA)) +
  theme_bw()

strengths_ttest <- strengths %>%
  group_by(net) %>%
  nest() %>%
  mutate(
    ttest = map(data, ~t.test(.x$strength2)),
    estimate = map_dbl(ttest, ~.x$estimate),
    p = map_dbl(ttest, ~.x$p.value)
  ) %>%
  add_column(p_adj = p.adjust(.$p, method = "fdr"))

ggplot(strengths, aes(x = net, y = strength2, fill = net)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(0, NA)) +
  theme_bw() +
  labs(x = "Node strength (within-network)",
       caption = "All mean values are >0, p<.001 (FDR)")

aov1_wnet <- aov(strength2 ~ net, data = strengths)

TukeyHSD(aov1_wnet)
