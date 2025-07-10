library(tidyverse)

# Setup ====

test_1class_ps <- read_rds("hemisphere_model_ps.rds")
nets_ordered <- read_rds("nets_ordered.rds")

nets_squares <- read_rds("nets_squares.rds") %>%
  mutate(
    start_order = start_order - 0.5,
    end_order = end_order + 0.5
  )

nets_order <- select(nets_ordered, roi, order)

source("00-enrichment_functions.R")

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
  )

net_plot <- ggplot(bootstrapped1, aes(x = net1, y = net2, fill = prop)) +
  geom_tile() +
  geom_text(aes(label = sig_label)) +
  viridis::scale_fill_viridis(limits = c(0, NA)) +
  scale_y_discrete(limits = rev) +
  coord_equal() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = NULL, title = "(C)")

# Create publication plot

(ld1_plot + pepper_plot) + plot_layout(axes = "collect")

# Sensitivity analysis ====

## Excluded networks ====

sens1 <- pepper %>%
  summarize_tibble_by_net(networks = nets_ordered, drop = FALSE) %>%
  mutate(
    prop = replace(prop, prop == 0, NA)
  ) %>%
  remove_tri()

sens1_bs <- bootstrap_pmatrix(pepper, nets = nets_ordered, drop = FALSE) %>%
  remove_tri()

ggplot(sens1_bs, aes(x = net1, y = net2, fill = sig)) +
  geom_tile() +
  geom_text(aes(label = sig)) +
  viridis::scale_fill_viridis() +
  scale_y_discrete(limits = rev) +
  coord_equal() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

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
