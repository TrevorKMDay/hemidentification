library(tidyverse)
library(patchwork)
library(viridis)

# Setup ====

setwd("~/MyDrive/Projects/hemisphere_fingerprinting/code/enrichment/")

source("00-enrichment_functions.R")

test_4class_ps <- read_rds("4class_model_ps.rds")
nets_ordered <- read_rds("nets_ordered.rds")

# Basic node visualization ====

high_value_p001 <- test_4class_ps %>%
  filter(
    p < .001
  ) %>%
  separate_wider_delim(feature, "_", names = c("ROI1", "ROI2"))

high_value_nodes_p001 <- high_value_p001 %>%
  select(name, ROI1, ROI2) %>%
  pivot_longer(c(ROI1, ROI2), names_to = "roi", values_to = "rois") %>%
  group_by(name, rois) %>%
  summarize(
    n = n()
  ) %>%
  arrange(name, desc(n))

# Plot all nodes that appear more than once in all three LDs
ggplot(filter(high_value_nodes_p001, n > 2), aes(y = rois, x = name)) +
  geom_tile(aes(fill = n)) +
  viridis::scale_fill_viridis() +
  theme_bw()

# Matrices ====

test_4class_wide <- test_4class_ps %>%
  pivot_wider(names_from = name, values_from = p, id_cols = feature) %>%
  separate_wider_delim(feature, "_", names = c("ROI1", "ROI2"))

# These are symmetric matrices of p values, in long format
g1 <- graph_to_tibble(test_4class_wide, "LD1")
g2 <- graph_to_tibble(test_4class_wide, "LD2")
g3 <- graph_to_tibble(test_4class_wide, "LD3")

p1 <- ggplot(g1, aes(x = ROI1, y = ROI2, fill = value < .05)) +
  geom_tile() +
  scale_x_discrete(limits = nets_ordered$roi, labels = NULL) +
  scale_y_discrete(limits = nets_ordered$roi, labels = NULL) +
  scale_fill_manual(values = c("white", "black"), na.value = "grey",) +
  coord_equal() +
  labs(title = "LD1")

p2 <- ggplot(g2, aes(x = ROI1, y = ROI2, fill = value < .05)) +
  geom_tile() +
  viridis::scale_fill_viridis() +
  scale_x_discrete(limits = nets_ordered$roi, labels = NULL) +
  scale_y_discrete(limits = nets_ordered$roi, labels = NULL) +
  scale_fill_manual(values = c("white", "black"), na.value = "grey") +
  labs(title = "LD2")

p3 <- ggplot(g3, aes(x = ROI1, y = ROI2, fill = value < .05)) +
  geom_tile() +
  scale_x_discrete(limits = nets_ordered$roi, labels = NULL) +
  scale_y_discrete(limits = nets_ordered$roi, labels = NULL) +
  scale_fill_manual(values = c("white", "black"), na.value = "grey") +
  labs(title = "LD3")

(p1 + p2 + p3) + plot_layout(axes = "collect", guides = "collect")

g4 <- g1 %>%
  rename(
    value1 = value
  ) %>%
  left_join(
    rename(g2, value2 = value)
  ) %>%
  left_join(
    rename(g3, value3 = value)
  ) %>%
  mutate(
    across(starts_with("value"), ~ .x < .05, .names = "{.col}_sig"),
    nsig = value1_sig + value2_sig + value3_sig
  )

ggplot(g4, aes(x = ROI1, y = ROI2, fill = factor(nsig))) +
  geom_tile() +
  scale_x_discrete(limits = nets_ordered$roi, labels = NULL) +
  scale_y_discrete(limits = nets_ordered$roi, labels = NULL) +
  scale_fill_manual(values = c("white", "grey", "orange", "red")) +
  coord_equal() +
  labs(title = "Frequency across all LDs")

# Results by network ====

bad_nets <- c("OAN", "pMN", "vMN")

g1_bynet <- g1 %>%
  mutate(
    sig = value < .05
  ) %>%
  summarize_tibble_by_net(., nets_ordered, drop = bad_nets) %>%
  remove_tri()

g2_bynet <- g2 %>%
  mutate(
    sig = value < .05
  ) %>%
  summarize_tibble_by_net(., nets_ordered, drop = bad_nets) %>%
  remove_tri()

g3_bynet <- g3 %>%
  mutate(
    sig = value < .05
  ) %>%
  summarize_tibble_by_net(., nets_ordered, drop = bad_nets) %>%
  remove_tri()

ggplot(g1_bynet, aes(x = net1, y = net2, fill = prop)) +
  geom_tile() +
  geom_text(aes(label = sig)) +
  scale_y_discrete(limits = rev) +
  scale_fill_viridis(limits = c(0, 0.15)) +
  coord_equal() +
  theme_bw()

ggplot(g2_bynet, aes(x = net1, y = net2, fill = prop)) +
  geom_tile() +
  geom_text(aes(label = sig)) +
  scale_y_discrete(limits = rev) +
  scale_fill_viridis(limits = c(0, 0.15)) +
  coord_equal() +
  theme_bw()

ggplot(g3_bynet, aes(x = net1, y = net2, fill = prop)) +
  geom_tile() +
  geom_text(aes(label = sig)) +
  scale_y_discrete(limits = rev) +
  scale_fill_viridis(limits = c(0, 0.15)) +
  coord_equal() +
  theme_bw()

g1ps <- g1 %>%
  mutate(
    sig = value < .05
  )

g2ps <- g2 %>%
  mutate(
    sig = value < .05
  )

g3ps <- g3%>%
  mutate(
    sig = value < .05
  )

g1_results <- bootstrap_pmatrix(g1ps, nets = nets_ordered) %>%
  remove_tri()

g2_results <- bootstrap_pmatrix(g2ps, nets = nets_ordered) %>%
  remove_tri()

g3_results <- bootstrap_pmatrix(g3ps, nets = nets_ordered) %>%
  remove_tri()

results <- bind_rows(g1_results, g2_results, g3_results, .id = "LD") %>%
  mutate(
    sig = replace(sig, sig == 0, NA),
    prop = replace(prop, prop == 0, NA),
    sig_label = case_when(
      pvalue < .001 ~ "***",
      pvalue < .01 ~ "**",
      pvalue < .05 ~ "*",
      TRUE ~ ""
    ),
  ) %>%
  filter(
    !(net1 %in% c("OAN", "pMN", "vMN")),
    !(net2 %in% c("OAN", "pMN", "vMN"))
  )

ggplot(results, aes(x = net1, y = net2, fill = prop)) +
  geom_tile() +
  geom_text(aes(label = sig_label)) +
  scale_fill_viridis(limits = c(0, NA)) +
  scale_y_discrete(limits = rev) +
  facet_wrap(vars(LD)) +
  coord_equal() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = NULL)
