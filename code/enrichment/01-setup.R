library(tidyverse)
library(reticulate)
library(qs)
library(infer)
library(igraph)
library(patchwork)
library(Matrix)

# Setup ====

setwd("~/MyDrive/Projects/hemisphere_fingerprinting/code/")
select <- dplyr::select

networks0 <- read_csv(file = paste0("../data/ColeAnticevicNetPartition/",
                                    "all_parcels_to_networks.csv"),
                      show_col_types = FALSE)

mismatched_parcels <- networks0 %>%
  filter(
    L != R
  )

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

nets_ordered <- nets %>%
  arrange(net_short, roi) %>%
  mutate(
    order = row_number(),
  ) %>%
  group_by(network) %>%
  mutate(
    order_in_net = row_number()
  )

write_rds(nets_ordered, "nets_ordered.rds")

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

second_model <- read_csv("modeling/lda_full/full_class_scalings.csv",
                         show_col_types = FALSE)

cor(true_model$LD1, second_model$LD1) # 0.86
cor(true_model$LD1, second_model$LD1, method = "s") # 0.85

scnd_model_top12 <- second_model %>%
  arrange(desc(LD1)) %>%
  mutate(
    n = row_number()
  ) %>%
  filter(
    n <= 12
  )

intersect(true_model_top12$feature, scnd_model_top12$feature) # 6 / 12


## Overlap ====

intersect_features <- function(model1, model2, n) {

  model1 <- arrange(model1, desc(LD1))
  model2 <- arrange(model2, desc(LD1))

  overlap <- intersect(model1$feature[1:n], model2$feature[1:n])

  return(length(overlap))

}

n_overlapping <- tibble(n = seq(10, nrow(true_model), by = 200)) %>%
  mutate(
    intersect = map_int(n, ~intersect_features(true_model, second_model, .x)),
    pct = intersect / n
  )

ggplot(n_overlapping, aes(x = n, y = pct)) +
  geom_point() +
  theme_bw() +
  labs(title = "# overlapping features in top n features")

## Read in bootstrapping -====

### 1-class model

bs_file <- "lda_full_bootstrapped.rds"

if (!file.exists(bs_file)) {

  bs_f <- list.files("modeling/lda_full/", "full_bs[0-9]*_scalings.csv",
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

  bs_fs2 <- list.files("modeling/lda_full/", "full_class_bs[0-9]*_scalings.csv",
                     full.names = TRUE)

  bs2_0 <- tibble(f = bs_fs2) %>%
    mutate(
      data = map(f, read_csv, show_col_types = FALSE, .progress = TRUE),
      n = row_number()
    )

  bs2 <- bs2_0 %>%
    select(n, data) %>%
    unnest(data)

  qsave(bs2, bs_file2)

} else {

  bs <- qread(bs_file, nthreads = 2)

}

# Run tests ====

## 1-class =====

test_1class <- bs %>%
  group_by(feature) %>%
  nest() %>%
  left_join(true_model, join_by(feature)) %>%
  mutate(

    p1 = map2_dbl(data, LD1,
                  ~get_p_value(tibble(LD1 = .x$LD1), .y,
                               direction = "two_sided")[[1]],
                 .progress = TRUE),

  )

test_1class_ps <- test_1class %>%
  select(feature, starts_with("p")) %>%
  pivot_longer(-feature, values_to = "p") %>%
  mutate(
    log_p = log10(p) %>%
      replace(., is.infinite(.), -4)
  )

## 1-class =====

scnd_model2 <- second_model %>%
  pivot_longer(-feature, values_to = "true_value")

test_4class <- bs2 %>%
  select(-n) %>%
  pivot_longer(-feature) %>%
  group_by(feature, name) %>%
  nest() %>%
  left_join(scnd_model2, join_by(feature, name)) %>%
  mutate(

    p1 = map2_dbl(data, true_value,
                  ~get_p_value(tibble(value = .x$value), .y,
                               direction = "two_sided")[[1]],
                  .progress = TRUE),

  )

test_4class_ps <- test_4class %>%
  select(feature, name, starts_with("p")) %>%
  mutate(
    log_p = log10(p1) %>%
      replace(., is.infinite(.), -4)
  )

write_rds(test_1class_ps, "hemisphere_model_ps.rds")
write_rds(test_4class_ps, "4class_model_ps.rds")

## Features p < .001 ====

high_value_p001 <- test_ps %>%
  filter(
    name == "p1",
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

length(unique(high_value_p001$rois))

## Features p < .01 ====

high_value_p01 <- test_ps %>%
  filter(
    name == "p1",
    p <= .01
  ) %>%
  separate_wider_delim(feature, "_", names = c("ROI1", "ROI2"))

n_hv_nodes_p01 <- length(unique(high_value_nodes_p01$rois))

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


ggplot(test_ps, aes(x = -1 * log_p, color = name)) +
  geom_density() +
  theme_bw()

lt_p05 <- sum(test2$p1 < .05)
lt_p01 <- sum(test2$p1 < .01)
lt_p001 <- sum(test2$p1 < .001)

# Convert to corr matrix ====



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

graph_to_matrix <- function(g_input, col) {

  g <- graph_from_data_frame(g_input[, 1:2], directed = FALSE)
  edge.attributes(g)$weight <- as.vector(g_input[col])

  gmatrix <- as_adjacency_matrix(g, attr = "weights", type = "both",
                                 sparse = FALSE) %>%
    as.matrix()

  return(gmatrix)

}

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

mismatched_parcels2 <- mismatched_parcels %>%
  left_join(nets_ordered, by = join_by(label2 == roi))

p2 +
  geom_vline(data = mismatched_parcels2, aes(xintercept = order),
             color = "red", alpha = 0.5) +
  geom_hline(data = mismatched_parcels2, aes(yintercept = order),
             color = "red", alpha = 0.5)

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
      sig = sum(sig, na.rm = TRUE),
      .groups = "keep"
    ) %>%
    mutate(
      prop = sig / n
    )

  return(result)

}

true_nets <- summarize_matrix_by_net(pmatrix, nets_ordered) %>%
  mutate(
    prop = replace(prop, prop == 0, NA)
  )

ggplot(true_nets, aes(x = net1, y = net2, fill = prop * 100)) +
  geom_tile() +
  coord_equal() +
  scale_y_discrete(limits = rev(nets_squares$net_short)) +
  viridis::scale_fill_viridis(limits = c(0, NA)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Network 1", y = "Network 2",
       caption = "Networks with no significant connections are greyed out")

## Create null distribution ====


shuffle_matrix <- function(mat) {

  total_n <- sum(mat, na.rm = TRUE) / 2

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

  # Replace names from original matrix
  colnames(mat_final) <- colnames(mat)
  rownames(mat_final) <- rownames(mat)

  return(mat_final)

}

bootstrap_pmatrix <- function(pmatrix, n = 10000) {

  message("Creating shuffled matrices")

  bootstrap1 <- tibble(rep = 1:10000) %>%
    mutate(

      mat = map(rep, ~shuffle_matrix(pmatrix), .progress = TRUE),

    )

  message("Summarizing reps ... ")

  bootstrap2 <- bootstrap1 %>%
    mutate(
      result = map(mat, ~summarize_matrix_by_net(.x, nets_ordered))
    ) %>%
    select(-mat) %>%
    unnest(result) %>%
    group_by(net1, net2, n) %>%
    nest() %>%
    left_join(true_nets, by = join_by(net1, net2, n))

  # Move p value calculation out to speed up reproc
  # bootstrap2p <- bootstrap2 %>%
  #   mutate(
  #
  #     pvalue = map2_dbl(data, sig,
  #                       ~get_p_value(tibble(sig = .x$sig), .y,
  #                                    direction = "greater")[[1]],
  #                       .progress = TRUE)
  #
  #   )

}

bootstrap2_meanvals <- bootstrap2 %>%
  filter(
    !(net1 %in% c("OAN", "pMN", "vMN")),
    !(net2 %in% c("OAN", "pMN", "vMN")),
  ) %>%
  select(-n, -sig, -prop) %>%
  unnest(data) %>%
  group_by(net1, net2) %>%
  summarize(
    n = n(),
    mean_nsig = mean(sig),
    sd_nsig = sd(sig)
  ) %>%
  mutate(
    se_nsig = sd_nsig / sqrt(n),
    net_net = paste(net1, net2, sep = "_"),
    within_net = net1 == net2
  )


bootstrap3 <- bootstrap2p %>%
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

true_nets3_duplicate <- true_nets3 %>%
  bind_rows(
    rename(true_nets3, net1 = net2, net2 = net1)
  )

ggplot(bootstrap2_meanvals, aes(x = net2)) +
  geom_pointrange(aes(y = mean_nsig,
                      ymin = mean_nsig - 1.96 * sd_nsig,
                      ymax = mean_nsig + 1.96 * sd_nsig,
                      color = within_net),
                  size = 0.25) +
  geom_point(
    data = true_nets3_duplicate,
    aes(y = n_sig),
    size = 1.5
  ) +
  geom_text(
    data = true_nets3_duplicate,
    aes(label = label),
    y = 40
  ) +
  scale_y_continuous(limits = c(NA, 42)) +
  facet_wrap(vars(net1), scales = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  labs(x = "Network 2", y = "# sig. features (1.96 × ±SD)")


png("plots/all_enrichment_plots.png", width = 6.5, height = 5, units = "in",
    res =300)

(p1 + p2) / (enhanced_plot) + plot_layout(heights = c(10, 10))

dev.off()
