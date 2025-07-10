library(tidyverse)
library(reticulate)

# Setup ====

select <- dplyr::select

setwd("~/MyDrive/Projects/hemisphere_fingerprinting/code/")

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

# Results ====

true_model <- read_csv("modeling/lda_full/full_scalings.csv",
                       show_col_types = FALSE)

true_model_top12 <- true_model %>%
  arrange(desc(LD1)) %>%
  mutate(
    n = row_number()
  ) %>%
  filter(
    n <= 12
  )

bs_f <- list.files("modeling/lda_full/", "full_bs.*_scalings.csv",
                   full.names = TRUE)

bs0 <- tibble(f = bs_f) %>%
  mutate(
    data = map(f, read_csv, show_col_types = FALSE, .progress = TRUE),
    n = row_number()
  )

bs <- bs0 %>%
  select(n, data) %>%
  unnest(data)


ggplot(test, aes(x = LD1)) +
  geom_density() +
  geom_vline(
    data = true_model_top12,
    aes(xintercept = LD1),
    color = "red"
  ) +
  facet_wrap(vars(feature)) +
  theme_bw()

library(infer)

test2 <- bs %>%
  group_by(feature) %>%
  nest() %>%
  left_join(true_model, join_by(feature)) %>%
  mutate(

    p1 = map2_dbl(data, LD1,
                  ~get_p_value(tibble(LD1 = .x$LD1), .y,
                               direction = "two_sided")[[1]],
                 .progress = TRUE),

    p2 = map2_dbl(data, LD1,
                  ~get_p_value(tibble(LD1 = abs(.x$LD1)), abs(.y),
                               direction = "greater")[[1]],
                  .progress = TRUE),

  )

test_ps <- test2 %>%
  select(feature, p1, p2) %>%
  pivot_longer(-feature, values_to = "p") %>%
  mutate(
    log_p = log10(p) %>%
      replace(., is.infinite(.), -4)
  )

high_value <- test_ps %>%
  filter(
    name == "p1",
    p < .001
  ) %>%
  separate_wider_delim(feature, "_", names = c("ROI1", "ROI2"))

tibble(rois = c(high_value$ROI1, high_value$ROI2)) %>%
  group_by(rois) %>%
  summarize(
    n = n()
  ) %>%
  arrange(desc(n))

ggplot(test_ps, aes(x = -1 * log_p, color = name)) +
  geom_density() +
  theme_bw()

lt_p05 <- sum(test2$p1 < .05)
lt_p01 <- sum(test2$p1 < .01)
lt_p001 <- sum(test2$p1 < .001)

test3 <- test2 %>%
  ungroup() %>%
  filter(
    p1 < .00001 | p2 < .00001
  )

ggplot(filter(bs, feature %in% test3$feature), aes(x = abs(LD1))) +
  geom_density() +
  geom_vline(
    data = test3,
    aes(xintercept = abs(LD1)),
    color = "red"
  ) +
  scale_x_continuous(breaks = c(0, 0.05, 0.1)) +
  facet_wrap(vars(feature)) +
  theme_bw()

ggplot(filter(bs, feature %in% test3$feature), aes(x = LD1)) +
  geom_density() +
  geom_vline(
    data = test3,
    aes(xintercept = LD1),
    color = "red"
  ) +
  scale_x_continuous(breaks = c(-0.1, 0, 0.1)) +
  facet_wrap(vars(feature)) +
  theme_bw()

# Convert to corr matrix ====

library(igraph)

nets_ordered <- nets %>%
  arrange(net_short, roi) %>%
  mutate(
    order = row_number(),
  ) %>%
  group_by(network) %>%
  mutate(
    order_in_net = row_number()
  )

nets_squares <- nets_ordered %>%
  filter(
    order_in_net == 1 | order_in_net == max(order_in_net)
  ) %>%
  mutate(
    startend = if_else(order_in_net == 1, "start", "end")
  ) %>%
  pivot_wider(id_cols = c(network, net_short), names_from = startend,
              values_from = roi) %>%
  mutate(
    start_order = match(start, nets_ordered$roi),
    end_order = match(end, nets_ordered$roi)
  )

g_input <- test2 %>%
  select(feature, LD1, p1) %>%
  separate_wider_delim(feature, "_", names = c("ROI1", "ROI2")) %>%
  mutate(
    across(c(ROI1, ROI2), ~factor(.x, levels = nets_ordered$roi))
  ) %>%
  arrange(ROI1, ROI2)

graph <- graph_from_data_frame(g_input, directed = FALSE)
E(graph)$weights <- g_input$LD1

gmatrix <- as_adjacency_matrix(graph, attr = "weights", type = "both") %>%
  as.matrix()

pgraph <- graph
E(pgraph)$weights <- g_input$p1 < .01

pmatrix <- as_adjacency_matrix(pgraph, attr = "weights", type = "both") %>%
  as.matrix()

p2 <- pmatrix %>%
  as_tibble(rownames = "ROI1") %>%
  pivot_longer(-ROI1, names_to = "ROI2", values_to = "sig")

# gmatrix[upper.tri(gmatrix, diag = TRUE)] <- NA

g2 <- gmatrix %>%
  as_tibble(rownames = "ROI1") %>%
  pivot_longer(-ROI1, names_to = "ROI2")  %>%
  left_join(p2) %>%
  mutate(
    value = replace(value, value == 0, NA),
    ROI1_order = match(ROI1, nets_ordered$roi),
    ROI2_order = match(ROI2, nets_ordered$roi)
  )

labels <- c("", nets_squares$net_short) %>%
  replace(., . %in% c("OAN", "pMN", "vMN"), "")

## Colored plot ====

p1 <- ggplot(NULL) +
  geom_tile(
    data = g2,
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
    aes(
        # These adjustments stop AudN from starting at 0, off the graph
        xmin = if_else(start_order - 1 < 1, 1, start_order - 1),
        xmax = end_order,
        ymin = if_else(start_order - 1 < 1, 1, start_order - 1),
        ymax = end_order,
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

p2 <- ggplot(NULL) +
  geom_tile(
    data = g2,
    aes(x = ROI1_order, y = ROI2_order, fill = sig, alpha = sig),
    linewidth = 0.2
  ) +
  scale_x_continuous(labels = labels,
                     breaks = c(0.5, nets_squares$end_order + 0.5),
                     minor_breaks = NULL) +
  scale_y_continuous(labels = NULL,
                     breaks = c(0.5, nets_squares$end_order + 0.5),
                     minor_breaks = NULL,
                     trans = "reverse") +
  scale_fill_manual(values = c("white", "black"), guide = "none") +
  geom_rect(
    data = nets_squares,
    aes(
      # These adjustments stop AudN from starting at 0, off the graph
      xmin = if_else(start_order - 1 < 1, 0.5, start_order - 0.5),
      xmax = end_order + 0.5,
      ymin = if_else(start_order - 1 < 1, 0.5, start_order - 0.5),
      ymax = end_order + 0.5,
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

p1 + p2 + plot_layout(axis_titles = "collect")

# Bootstrap models ====

summarize_matrix_by_net <- function(mat, networks) {

  long <- mat %>%
    as_tibble(rownames = "ROI1") %>%
    pivot_longer(-ROI1, names_to = "ROI2", values_to = "sig") %>%
    mutate(
      net1 = networks$net_short[match(ROI1, networks$roi)],
      net2 = networks$net_short[match(ROI2, networks$roi)],
    )

  result <- long %>%
    group_by(net1, net2) %>%
    summarize(
      n = n(),
      sig = sum(sig),
      .groups = "keep"
    ) %>%
    mutate(
      p = sig / n
    )

  return(result)

}

true_nets <- summarize_matrix_by_net(pmatrix, nets_ordered) %>%
  mutate(
    p = replace(p, p == 0, NA)
  )

ggplot(true_nets, aes(x = net1, y = net2, fill = p)) +
  geom_tile() +
  coord_equal() +
  scale_y_discrete(limits = rev(nets_squares$net_short)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Network 1", y = "Network 2",
       caption = "Networks with no significant connections are greyed out")

## Create null distribution ====

library(Matrix)

shuffle_matrix <- function(mat) {

  total_n <- sum(mat) / 2

  mat_lower <- mat
  mat_lower[upper.tri(mat_lower, diag = TRUE)] <- NA

  v <- as.vector(mat_lower)

  v2 <- rep(NA, length(v))
  v2[!is.na(v)] <- FALSE

  new_random_lower <- sample(which(!is.na(v)), size = total_n)
  v2[new_random_lower] <- TRUE

  mat_final <- matrix(v2, nrow = nrow(mat_lower), ncol = ncol(mat_lower))

  mat_final[upper.tri(mat_final)] <- t(mat_final)[upper.tri(mat_final)]
  mat_final[diag(mat_final)] <- FALSE

  return(mat_final)

}

bootstrap1 <- tibble(rep = 1:10000) %>%
  mutate(
    mat = map(rep, ~shuffle_matrix(pmatrix),
              .progress = TRUE)
  )

bootstrap2 <- bootstrap1 %>%
  mutate(
    result = map(mat, ~summarize_matrix_by_net(.x, nets_ordered))
  ) %>%
  select(-mat) %>%
  unnest(result) %>%
  group_by(net1, net2, n) %>%
  nest() %>%
  left_join(true_nets) %>%
  mutate(
    pvalue = map2_dbl(data, sig,
                      ~get_p_value(tibble(sig = .x$sig), .y,
                                   direction = "greater")[[1]],
                      .progress = TRUE)
  )

bootstrap3 <- bootstrap2 %>%
  select(net1, net2, pvalue) %>%
  arrange(pvalue) %>%
  mutate(
    label = case_when(
      pvalue < .001 ~ "***",
      pvalue < .01 ~ "**",
      pvalue < .05 ~ "*",
      TRUE ~ ""
    )
  ) %>%
  filter(
    !(net1 %in% c("OAN", "pMN", "vMN")),
    !(net2 %in% c("OAN", "pMN", "vMN"))
  )

# Turn all values into matrix to remove upper triangle, then put back in long
# format

true_nets2 <- true_nets %>%
  filter(
    !(net1 %in% c("OAN", "pMN", "vMN")),
    !(net2 %in% c("OAN", "pMN", "vMN"))
  ) %>%
  pivot_wider(id_cols = net1, names_from = net2, values_from = sig) %>%
  column_to_rownames("net1")

true_nets2[lower.tri(true_nets2)] <- NA

true_nets3 <- true_nets2 %>%
  as_tibble(rownames = "net1") %>%
  pivot_longer(-net1, names_to = "net2", values_to = "n_sig") %>%
  na.omit() %>%
  left_join(select(true_nets, net1, net2, n)) %>%
  mutate(
    prop = if_else(n_sig / n > 0, 100 * n_sig / n, NA)
  ) %>%
  left_join(bootstrap3)

enhanced_plot <- ggplot(true_nets3, aes(x = net1, y = net2)) +
  geom_tile(aes(fill = prop)) +
  geom_text(aes(label = label)) +
  coord_equal() +
  scale_y_discrete(limits = rev) +
  viridis::scale_fill_viridis(limits = c(0, NA), na.value = "grey75",
                              labels = paste0(0:4, "%")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Network 1", y = "Network 2", title = "(C)")

png("plots/all_enrichment_plots.png", width = 6.5, height = 5, units = "in",
    res =300)

(p1 + p2) / (enhanced_plot) + plot_layout(heights = c(10, 10))

dev.off()
