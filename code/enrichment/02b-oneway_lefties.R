library(tidyverse)
library(patchwork)
library(viridis)
library(qs)

# Setup ====

setwd("~/MyDrive/Projects/hemisphere_fingerprinting/code/enrichment/")

source("00-enrichment_functions.R")

test_1class_ps <- read_rds("hemisphere_lefies_ps.rds")
nets_ordered <- read_rds("nets_ordered.rds")

nets_squares <- read_rds("nets_squares.rds") %>%
  mutate(
    start_order = start_order - 0.5,
    end_order = end_order + 0.5
  )

nets_order <- select(nets_ordered, roi, order)

## Features p < .001 ====

high_value_p001 <- test_1class_ps %>%
  filter(
    p < .001
  ) %>%
  separate_wider_delim(feature, "_", names = c("ROI1", "ROI2"))

high_value_nodes_p001 <- tibble(rois = c(high_value_p001$ROI1,
                                         high_value_p001$ROI2)) %>%
  group_by(rois) %>%
  summarize(
    n = n()
  ) %>%
  arrange(desc(n))

length(unique(high_value_nodes_p001$rois))

## Features p < .01 ====

high_value_p01 <- test_1class_ps %>%
  filter(
    p <= .01
  ) %>%
  separate_wider_delim(feature, "_", names = c("ROI1", "ROI2"))

n_hv_nodes_p01 <- length(unique(high_value_nodes_p001$rois))

high_value_nodes_p01 <- tibble(rois = c(high_value_p01$ROI1,
                                        high_value_p01$ROI2)) %>%
  group_by(rois) %>%
  summarize(
    n = n()
  ) %>%
  mutate(
    p = round(n / n_hv_nodes_p01, 3)
  ) %>%
  arrange(desc(n))


ggplot(test_1class_ps, aes(x = -1 * log_p)) +
  geom_density() +
  geom_vline(xintercept = -1 * log(c(0.05, 0.01, 0.001))) +
  theme_bw()

lt_p05 <- sum(test_1class_ps$p < .05)
lt_p01 <- sum(test_1class_ps$p < .01)
lt_p001 <- sum(test_1class_ps$p < .001)

# Visualizations ====

# raw p-value ====

g_input <- test_1class_ps %>%
  select(feature, LD1, p) %>%
  separate_wider_delim(feature, "_", names = c("ROI1", "ROI2")) %>%
  mutate(
    across(c(ROI1, ROI2), ~factor(.x, levels = nets_ordered$roi))
  ) %>%
  arrange(ROI1, ROI2)

g1 <- graph_to_tibble(g_input, "LD1") %>%
  mutate(
    ROI1_order = map_int(ROI1, ~ nets_order$order[nets_order$roi == .x]),
    ROI2_order = map_int(ROI2, ~ nets_order$order[nets_order$roi == .x])
  )

# pepper plot ====

p_input <- g_input %>%
  mutate(
    sig = p < .01
  )

pepper <- graph_to_tibble(p_input, "sig") %>%
  mutate(
    ROI1_order = map_int(ROI1, ~ nets_order$order[nets_order$roi == .x]),
    ROI2_order = map_int(ROI2, ~ nets_order$order[nets_order$roi == .x])
  ) %>%
  rename(
    sig = value
  )

## Colored plot ====

labels <- c("", nets_squares$net_short[1:6], "", nets_squares$net_short[7:8],
            "", "")

ld1_plot <- ggplot(NULL) +
  geom_tile(
    data = g1,
    aes(x = ROI1_order, y = ROI2_order, fill = value),
    linewidth = 0.2
  ) +
  scale_x_continuous(labels = labels,
                     breaks = c(1, nets_squares$end_order),
                     minor_breaks = NULL) +
  scale_y_continuous(labels = labels,
                     breaks = c(1, nets_squares$end_order),
                     minor_breaks = NULL,
                     trans = "reverse") +
  scale_fill_gradient2(na.value = NA, limits = c(-0.3, 0.3)) +
  geom_rect(
    data = nets_squares,
    aes(xmin = start_order, xmax = end_order,
        ymin = start_order, ymax = end_order,
        color = net_short),
    fill = NA
  ) +
  coord_equal() +
  theme_bw() +
  theme(axis.ticks = NULL, axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  labs(x = "ROI 1", y = "ROI 2", color = "Network", fill = "Scaling",
       title = "(A)") +
  guides(color = guide_legend(ncol = 2))

## Pepper plot =====

pepper_plot <- ggplot(NULL) +
  geom_tile(
    data = pepper,
    aes(x = ROI1_order, y = ROI2_order, fill = sig, alpha = sig),
    linewidth = 0.2
  ) +
  scale_x_continuous(labels = labels,
                     breaks = c(1, nets_squares$end_order),
                     minor_breaks = NULL) +
  scale_y_continuous(labels = labels,
                     breaks = c(1, nets_squares$end_order),
                     minor_breaks = NULL,
                     trans = "reverse") +
  scale_fill_manual(values = c("white", "black"), guide = "none",
                    na.value = "white") +
  geom_rect(
    data = nets_squares,
    aes(xmin = start_order, xmax = end_order,
        ymin = start_order, ymax = end_order,
        color = net_short),
    fill = NA
  ) +
  coord_equal() +
  theme_bw() +
  theme(axis.ticks = NULL, axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "ROI 1", y = NULL, color = "Network", fill = "Sig.",
       alpha = "Sig.",
       title = "(B)")+
  guides(color = guide_legend(ncol = 2), alpha = "none")

ld1_plot + pepper_plot + plot_layout(axis_titles = "collect")

# Which connections matter ====

sa_rank <- read_csv("https://raw.githubusercontent.com/PennLINC/S-A_ArchetypalAxis/refs/heads/main/Glasser360_MMP/Sensorimotor_Association_Axis_AverageRanks.csv") %>%
  select(-brains.average.rank)

hc <- qread("../modeling/inputs/hemiconnectome.rds")

connection_strengths <- hc %>%
  select(sub, hemi, contains("_"), -age_group) %>%
  pivot_longer(-c(sub, hemi)) %>%
  pivot_wider(names_from = hemi) %>%
  group_by(name) %>%
  summarize(
    mean_RH = mean(RH),
    mean_LH = mean(LH),
    mean_diff = mean(LH - RH)
  ) %>%
  separate_wider_delim(name, delim = "_", names = c("ROI1", "ROI2"))

important_connections <- pepper %>%
  filter(
    sig
  ) %>%
  left_join(nets_order, by = join_by(ROI1 == roi)) %>%
  left_join(nets_order, by = join_by(ROI2 == roi), suffix = c("1", "2")) %>%
  select(-contains("order")) %>%
  left_join(connection_strengths) %>%
  filter(
    !is.na(mean_diff)
  ) %>%
  mutate(
    conn = paste(ROI1, ROI2, sep = "_")
  ) %>%
  left_join(sa_rank, by = join_by(ROI1 == region)) %>%
  left_join(sa_rank, by = join_by(ROI2 == region), suffix = c("1", "2"))

ggplot(important_connections, aes(x = final.rank1, y = mean_diff)) +
  geom_point() +
  scale_x_continuous(breaks = important_connections$final.rank1,
                     labels = important_connections$ROI1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


important_connections2 <- important_connections %>%
  filter(
    network1 %in% c("Default", "Frontoparietal", "Language") &
      network2 %in% c("Default", "Frontoparietal", "Language")
  )

table(c(important_connections2$ROI1, important_connections2$ROI2)) %>%
  sort()

ggplot(important_connections2, aes(x = conn, y = mean_diff)) +
  geom_point() +
  geom_hline(yintercept = 0, color = "red") +
  facet_grid(rows = vars(network1), cols = vars(network2),
             scales = "free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Bootstrap models ====

true_nets <- summarize_tibble_by_net(pepper, nets_ordered) %>%
  mutate(
    prop = replace(prop, prop == 0, NA)
  ) %>%
  remove_tri()

bootstrapped <- bootstrap_pmatrix(pepper, nets_ordered)

bootstrapped1 <- remove_tri(bootstrapped) %>%
  mutate(

    sig = replace(sig, sig == 0, NA),
    prop = replace(prop, prop == 0, NA),

    sig_label = case_when(
      pvalue < .001 ~ "***",
      pvalue < .01 ~ "**",
      pvalue < .05 ~ "*",
      TRUE ~ ""
    ),

    pvalue = replace(pvalue, pvalue == 0, 0.0003)

  )

bootstrapped1$p_fdr <- p.adjust(bootstrapped1$pvalue, method = "fdr")

bootstrapped1 %>%
  filter(
    pvalue < .05
  )

net_plot <- ggplot(bootstrapped1, aes(x = net1, y = net2, fill = prop)) +
  geom_tile() +
  geom_text(aes(label = sig_label)) +
  viridis::scale_fill_viridis(limits = c(0, NA)) +
  scale_y_discrete(limits = rev) +
  coord_equal() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = NULL, title = "Enrichment analysis, sinistrals only")

# Create publication plot

(ld1_plot + pepper_plot) + plot_layout(axes = "collect")

net_plot

# Sensitivity analysis ====

## Excluded networks ====

sens1 <- pepper %>%
  summarize_tibble_by_net(networks = nets_ordered, drop = FALSE) %>%
  mutate(
    prop = replace(prop, prop == 0, NA)
  ) %>%
  remove_tri()

sens1_bs <- bootstrap_pmatrix(pepper, nets = nets_ordered, drop = FALSE) %>%
  remove_tri() %>%
  mutate(

    sig = replace(sig, sig == 0, NA),

    sig_label = case_when(
      pvalue < .001 ~ "***",
      pvalue < .01 ~ "**",
      pvalue < .05 ~ "*",
      TRUE ~ ""
    ),
  )

ggplot(sens1_bs, aes(x = net1, y = net2, fill = sig)) +
  geom_tile() +
  geom_text(aes(label = sig_label)) +
  scale_fill_viridis(limits = c(0, NA)) +
  scale_y_discrete(limits = rev) +
  coord_equal() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Sensitivity: All networks")

## Flip-flopped ROIs ====

networks0 <- read_csv(file = paste0("../../data/ColeAnticevicNetPartition/",
                                    "all_parcels_to_networks.csv"),
                      show_col_types = FALSE)

mismatched_parcels <- networks0 %>%
  filter(
    L != R
  )

sens2 <- pepper %>%
  filter(
    !(ROI1 %in% mismatched_parcels$label2),
    !(ROI2 %in% mismatched_parcels$label2)
  )

sens2_bs <- bootstrap_pmatrix(sens2, nets = nets_ordered) %>%
  remove_tri()

sens2_bs1 <- sens2_bs %>%
  mutate(

    sig = replace(sig, sig == 0, NA),
    prop = replace(prop, prop == 0, NA),

    sig_label = case_when(
      pvalue < .001 ~ "***",
      pvalue < .01 ~ "**",
      pvalue < .05 ~ "*",
      TRUE ~ ""
    ),

    pvalue = replace(pvalue, pvalue == 0 , .0003)

  )

sens2_bs1$p_fdr <- p.adjust(sens2_bs1$pvalue, method = "fdr")

sens2_bs1 %>%
  filter(
    pvalue < .05
  )

ggplot(sens2_bs1, aes(x = net1, y = net2, fill = prop)) +
  geom_tile() +
  geom_text(aes(label = sig_label)) +
  viridis::scale_fill_viridis() +
  scale_y_discrete(limits = rev) +
  coord_equal() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Sensitivity: Mismatched ROIs removed")

