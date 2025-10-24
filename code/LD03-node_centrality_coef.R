library(tidyverse)
library(igraph)
library(corrr)
library(GGally)
library(ciftiTools)

# Setup ====

setwd("~/MyDrive/Projects/hemisphere_fingerprinting/code/")
ciftiTools.setOption('wb_path',
                     "/Users/tkmd/bin/workbench/bin_macosx64/wb_command/")

# Load data ====



hcp_cogs <- "https://bitbucket.org/dpat/tools/raw/master/REF/ATLASES/HCP-MMP1_UniqueRegionList.csv"

cogs <- read_csv(hcp_cogs, show_col_types = FALSE) %>%
  filter(
    LR == "L"
  ) %>%
  select(regionName, ends_with("cog")) %>%
  mutate(
    roi = str_remove(regionName, "_L") %>%
      str_replace("7Pl", "7PL")
  ) %>%
  select_all(~str_remove(., "-cog")) %>%
  select(roi, x, y, z)

dist <- expand_grid(cogs, select_all(cogs, ~paste0(., "2"))) %>%
  filter(
    roi != roi2
  ) %>%
  mutate(
    feature = paste0(roi, "_", roi2)
  ) %>%
  select(-roi, -roi2) %>%
  pivot_longer(-feature) %>%
  mutate(
    pair = if_else(str_detect(name, "2"), "1", "2"),
    name = str_remove(name, "2$")
  ) %>%
  pivot_wider(names_from = pair, names_prefix = "roi") %>%
  group_by(feature) %>%
  nest() %>%
  ungroup() %>%
  mutate(
    dist_vox = map_dbl(data, ~sqrt(sum((.x$roi1 - .x$roi2)^2))),
    dist_vox_Z = scale(dist_vox)[, 1]
  ) %>%
  select(-data)

ggplot(dist, aes(dist_vox)) +
  geom_density()

coef_files <- list.files("modeling/results/base", "*_coefs.csv",
                            full.names = TRUE)

coefs <- tibble(f = coef_files) %>%
  mutate(
    group = str_extract(f, "test-.") %>%
      str_remove("test-"),
    data = map(f, read_csv, show_col_types = FALSE)
  ) %>%
  unnest(data) %>%
  rename(
    feature = `...1`
  ) %>%
  select(-f)

write_rds(coefs, "coefs.rds")

to_pscalar <- read_table("degree_centrality_to_pscalar.txt", col_names = "value")

coefs_wide <- coefs %>%
  pivot_wider(id_cols = feature, names_from = group,
              values_from = ends_with("H"))

write_tsv(coefs_wide, "coefs_wide.txt")

# coefs_by_group <- coefs %>%
#   pivot_wider(id_cols = feature, names_from = group,
#               values_from = ends_with("H"))
#
# correlate(coefs_by_group) %>%
#   rplot()

coefs2 <- coefs %>%
  group_by(group) %>%
  mutate(
    # Don't center this bc 0 should be 0
    across(ends_with("H"), ~scale(.x, center = FALSE)[,1], .names = "{.col}_Z"),
    righty_Z = rightyRH_Z - rightyLH_Z,
    lefty_Z = leftyRH_Z - leftyLH_Z
  )

correlate(coefs2)

prominence <- function(v) {

  result <- rep(NA, length(v))
  for (i in 1:length(v)) {

    p <- v[i]
    rest <- v[-i]

    p2 <- p - mean(rest)

    result[i] <- p2

  }

  return(result)

}

coefs_prominence <- coefs %>%
  pivot_longer(-c(group, feature)) %>%
  group_by(feature) %>%
  mutate(
    Z = scale(value)[, 1],
    p = prominence(value)
  ) %>%
  arrange(feature, name, group)

ggplot(coefs_prominence, aes(x = abs(Z))) +
  geom_density(aes(color = name), alpha = 0.25) +
  theme_bw()

coefs_prom_wide <- coefs_prominence %>%
  pivot_wider(values_from = Z, id_cols = c(group, feature))

coefs_prominence %>%
  mutate(
    Z = abs(Z)
  ) %>%
  pivot_wider(values_from = Z, names_from = c(group, name),
              id_cols = feature) %>%
  correlate() %>%
  clipr::write_clip()

correlate(coefs_prom_wide)

coefs_mean <- coefs2 %>%
  select(feature, group, ends_with("_Z")) %>%
  pivot_longer(-c(feature, group)) %>%
  group_by(feature, name) %>%
  summarize(
    mean_zcoef = mean((value))
  ) %>%
  left_join(dist, join_by(feature))

# coefs_mean %>%
#   filter(
#     is.na(dist_vox)
#   )

coefs_mean %>%
  pivot_wider(values_from = mean_zcoef) %>%
  correlate()

ggplot(coefs_mean, aes(x = dist_vox, y = abs(mean_zcoef))) +
  geom_point(alpha = 0.05) +
  geom_smooth(method = "lm") +
  facet_wrap(vars(name))

coefs_distance_models <- coefs_mean %>%
  group_by(name) %>%
  nest() %>%
  mutate(
    lm = map(data, ~lm(mean_zcoef ~ dist_vox_Z, data = .x))
  )

## Networks ====

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

net_sizes <- nets %>%
  group_by(net_short) %>%
  summarize(
    n = n()
  )

# Descriptive plots ====

ggplot(coefs_mean, aes(mean_zcoef, color = name)) +
  geom_density()

coefs_mean_wide <- pivot_wider(coefs_mean, values_from = mean_zcoef)

# Coefs against all coefs
ggpairs(select(ungroup(coefs_mean_wide), -feature))

# Reduce by node ====

coefs_mean_rois <- coefs_mean %>%
  separate_wider_delim(feature, "_", names = c("ROI1", "ROI2")) %>%
  left_join(
    select(nets, roi, net_short), join_by(ROI1 == roi)
  ) %>%
  left_join(
    select(nets, roi, net_short), join_by(ROI2 == roi), suffix = c("1", "2")
  ) %>%
  mutate(
    name = str_remove(name, "_Z"),
    withinnet = (net_short1 == net_short2)
  )

ggplot(coefs_mean_rois, aes(x = name, y = mean_zcoef, fill = withinnet)) +
  geom_boxplot()

create_graph <- function(connections) {

  graph <- graph_from_data_frame(connections, directed = FALSE) %>%
    simplify()

  nodes <- names(V(graph))
  E(graph)$weight <- abs(connections$scaling)

  return(graph)

}

strength_from_graph <- function(graph) {

  nodes <- names(V(graph))
  eigen <- eigen_centrality(graph)

  strength_df <- tibble(node = nodes, strength = strength(graph),
                        eigen = eigen$vector) %>%
    arrange(node) %>%
    mutate(
      z = scale(strength)[, 1],
      q = percent_rank(strength)
    )

  return(strength_df)

}

coef_mean_rois2 <- coefs_mean_rois %>%
  group_by(name) %>%
  mutate(
    scaling = mean_zcoef
  ) %>%
  nest() %>%
  mutate(
    G = map(data, create_graph),
    str = map(G, strength_from_graph),

    withinnet = map(data, ~filter(.x, withinnet)),
    Gw = map(withinnet, create_graph),
    strw = map(Gw, strength_from_graph),
  )

# Extract strengths
#   Filter out networks we don't care about, and adjust strength by network size
strengths <- coef_mean_rois2 %>%
  select(name, str, strw) %>%
  unnest(c(str, strw), names_sep = "_") %>%
  select(-strw_node) %>%
  left_join(
    select(nets, roi, net_short),
    join_by(str_node == roi)
  ) %>%
  filter(
    !(net_short %in% c("OAN", "pMN", "vMN"))
  ) %>%
  left_join(net_sizes, by = join_by(net_short)) %>%
  mutate(
    strw_adj = strw_strength / n,
    strw_eigen_adj = strw_eigen / n
  ) %>%
  select(str_node, net_short, name, everything())

strengths_long <- strengths %>%
  select(name, str_node, net_short, ends_with("strength"), ends_with("adj"),
         ends_with("eigen")) %>%
  pivot_longer(c(ends_with("strength"), ends_with("adj"), ends_with("eigen")),
               names_to = "model",
               values_to = "strength") %>%
  mutate(
    model = case_match(
      model,
      "str_strength" ~ "all",
      "strw_strength" ~ "wnet",
      "strw_adj" ~ "wnet_adj",
      "str_eigen" ~ "all_eigen",
      "strw_eigen" ~ "wnet_eigen",
      "strw_eigen_adj" ~ "wnet_eigen_adj"
    ),
    diff_score = if_else(name %in% c("lefty", "righty"), "diff_score", "score")
  )

ggplot(strengths_long, aes(x = net_short, y = strength)) +
  geom_boxplot(aes(fill = name)) +
  facet_grid(rows = vars(model), cols = vars(diff_score), scales = "free") +
  theme_bw()

ggplot(strengths, aes(x = str_strength, y = strw_strength)) +
  geom_point(alpha = 0.25) +
  geom_smooth(method = "lm") +
  facet_wrap(vars(name)) +
  labs(x = "Strength (all)", y = "Strength (within-net)") +
  theme_bw()

ggplot(strengths, aes(x = str_strength, y = strw_adj)) +
  geom_point(alpha = 0.25) +
  geom_smooth(method = "lm") +
  facet_wrap(vars(name)) +
  labs(x = "Strength (all)", y = "Strength (within-net)") +
  theme_bw()

aov(strw_strength ~ net_short*name, data = strengths) %>%
  summary()

ggplot(strengths, aes(x = str_node, y = name, fill = str_strength)) +
  geom_tile() +
  viridis::scale_fill_viridis() +
  facet_wrap(vars(net_short), scales = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "ROI", y = "Hemi",
       title = "All strength")

ggplot(strengths, aes(x = str_node, y = name, fill = strw_strength)) +
  geom_tile() +
  viridis::scale_fill_viridis() +
  facet_wrap(vars(net_short), scales = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "ROI", y = "Hemi",
       title = "Within-net strength")

ggplot(strengths, aes(x = str_node, y = name, fill = strw_adj)) +
  geom_tile() +
  viridis::scale_fill_viridis() +
  facet_wrap(vars(net_short), scales = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "ROI", y = "Hemi",
       title = "Within-net strength, adj.")

ggplot(strengths, aes(x = str_node, y = name, fill = strw_eigen)) +
  geom_tile() +
  viridis::scale_fill_viridis() +
  facet_wrap(vars(net_short), scales = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "ROI", y = "Hemi",
       title = "Within-net strength, adj.")

strengths <- strengths %>%
  mutate(
    hand = str_remove(name, ".H$"),
    hemi = str_extract(name, ".H$")
  )

lm(strw_strength ~ net_short*hand*hemi, data = strengths) %>%
  summary()

lm(strw_adj ~ net_short*hand*hemi, data = strengths) %>%
  summary()

# Correlation between the four classes is .95 - .998
strengths_by_node <- strengths %>%
  pivot_wider(id_cols = c(str_node, net_short),
              names_from = name,
              values_from = c(str_strength, strw_strength, strw_adj,
                              str_eigen)) %>%
  mutate(
    mean_strength = select(., matches("str_strength_.*H")) %>%
      rowMeans(),
    mean_wstrength = select(., matches("strw_strength_.*H")) %>%
      rowMeans(),
    mean_strwadj = select(., matches("strw_adj_.*H")) %>%
      rowMeans(),
    mean_streigen = select(., matches("str_eigen_.*H")) %>%
      rowMeans(),
  ) %>%
  fastDummies::dummy_cols("net_short", omit_colname_prefix = TRUE)

ggplot(strengths_by_node, aes(x = net_short, y = mean_strength)) +
  geom_violin() +
  geom_jitter(height = 0, width = 0.1, alpha = 0.75) +
  theme_bw()

ggplot(strengths_by_node, aes(x = net_short, y = strw_strength_righty)) +
  geom_violin() +
  geom_jitter(height = 0, width = 0.1, alpha = 0.75) +
  theme_bw()

ggplot(strengths_by_node, aes(x = net_short, y = mean_wstrength)) +
  geom_violin() +
  geom_jitter(height = 0, width = 0.1, alpha = 0.75) +
  theme_bw()

ggplot(strengths_by_node, aes(x = net_short, y = mean_strwadj)) +
  geom_violin() +
  geom_jitter(height = 0, width = 0.1, alpha = 0.75) +
  theme_bw()

ggplot(strengths_by_node, aes(x = net_short, y = mean_streigen)) +
  geom_violin() +
  geom_jitter(height = 0, width = 0.1, alpha = 0.75) +
  theme_bw()


lm(mean_strength ~ AudN + CON + DAN + DMN + FPN + LN + SMN,
   data = strengths_by_node) %>%
  summary()

lm(mean_streigen ~ AudN + CON + DAN + DMN + FPN + LN + SMN,
   data = strengths_by_node) %>%
  summary()
