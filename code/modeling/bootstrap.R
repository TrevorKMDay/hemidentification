library(tidyverse)
library(mltools)

setwd("~/Google Drive/My Drive/Projects/hemisphere_fingerprinting/code/modeling")

# Networks ====

sizes <- read_rds("connection_sizes.rds") %>%
  rename(
    size = n
  )

files_networks <- list.files("results/networks/", pattern = "*network.*.csv",
                             full.names = TRUE)

results <- tibble(f = files_networks) %>%
  mutate(
    data = map(f, read_csv, show_col_types = FALSE,
               name_repair = "unique_quiet", .progress = TRUE),
    f = basename(f) %>%
      str_remove(".csv")
  )

results_big <- results %>%
  mutate(

    n = map_int(data, nrow),

    method = str_extract(f, "method-[^_]*"),
    outcome = str_extract(f, "outcome-[^_]*"),
    group = str_extract(f, "test-[^_]*"),
    input = str_extract(f, "data-.*"),

    hands = replace_na(str_extract(f, "hands-[^_]*"), "all"),
    hemis = replace_na(str_extract(f, "hemi-[^_]*"), "all"),

    across(c(outcome, method, group, hands, hemis, input),
           ~str_remove(.x, "^.*-"))

  ) %>%
  unnest(data) %>%
  mutate(
    ground_hand = str_remove(ground, ".H$"),
    ground_hemi = str_remove(predicted, "^(left|right)y"),
    pred_hand = str_remove(predicted, ".H$"),
    pred_hemi = str_remove(predicted, "^(left|right)y"),
    across(c(starts_with("ground"), starts_with("pred")), as.factor),
  ) %>%
  select(-`...1`, -f) %>%
  group_by(outcome, method, group, hands, hemis, input, n)

networks <- results_big %>%
  filter(
    !(input %in% c("base", "osampl", "trans")),
    !str_detect(input, "_")
  ) %>%
  nest() %>%
  select(outcome, method, group, hands, hemis, input, data, n) %>%
  mutate(

    confusion = map(data, ~as.data.frame.matrix(table(.x$ground, .x$predicted))),
    mcc = map_dbl(confusion, ~mcc(confusionM = .x)),

    confusion_hemi = map(data, ~as.data.frame.matrix(table(.x$ground_hemi,
                                                           .x$pred_hemi))),
    mcc_hemi = map_dbl(confusion_hemi, ~mcc(confusionM = .x)),

    confusion_hand = map(data, ~as.data.frame.matrix(table(.x$ground_hand,
                                                           .x$pred_hand))),
    mcc_hand = map_dbl(confusion_hand, ~mcc(confusionM = .x)),

  ) %>%
  mutate(
    net1 = input,
    net2 = input
  ) %>%
  left_join(sizes)

nets <- unique(networks$input)

png("plots/networks1.png", width = 5.5, height = 4, units = "in",
    res = 300)

ggplot(networks, aes(x = input, y = mcc, color = group)) +
  geom_point() +
  scale_y_continuous(limits = c(0, 1)) +
  facet_grid(cols = vars(method), rows = vars(hands)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Network", y = "MCC",
        title = "Four-way prediction trained on handedness groups")

dev.off()

ggplot(networks, aes(x = size, y = mcc)) +
  geom_point(aes(color = net1)) +
  geom_smooth() +
  scale_x_log10() +
  facet_grid(cols = vars(method), rows = vars(hands)) +
  theme_bw()

# ggplot(networks, aes(x = input, y = mcc_hemi, color = group)) +
#   geom_point() +
#   scale_y_continuous(limits = c(0, 1)) +
#   facet_grid(cols = vars(method), rows = vars(hands)) +
#   theme_bw() +
#   theme(axis.text = element_text(angle = 45, hjust = 1))

png("plots/networks_hands.png", width = 5.5, height = 4, units = "in",
    res = 300)

ggplot(filter(networks, hands == "all"),
       aes(x = input, y = mcc_hand, color = group)) +
  geom_point() +
  scale_y_continuous(limits = c(-0.2, 0.2)) +
  facet_grid(cols = vars(method)) +
  theme_bw() +
  theme(axis.text = element_text(angle = 45, hjust = 1))

dev.off()

two_networks <- results_big %>%
  filter(
    str_detect(input, "_")
  ) %>%
  nest() %>%
  select(outcome, method, group, hands, hemis, input, data, n) %>%
  separate_wider_delim(input, delim = "_", names = c("net1", "net2")) %>%
  mutate(

    confusion = map(data,
                    ~as.data.frame.matrix(table(.x$ground, .x$predicted))),
    mcc = map_dbl(confusion, ~mcc(confusionM = .x)),

    confusion_hemi = map(data, ~as.data.frame.matrix(table(.x$ground_hemi,
                                                           .x$pred_hemi))),
    mcc_hemi = map_dbl(confusion_hemi, ~mcc(confusionM = .x)),

    confusion_hand = map(data, ~as.data.frame.matrix(table(.x$ground_hand,
                                                           .x$pred_hand))),
    mcc_hand = map_dbl(confusion_hand, ~mcc(confusionM = .x)),

  ) %>%
  left_join(sizes)

png("plots/between_networks.png", width = 5.5, height = 4, units = "in",
    res = 300)

two_networks_summary <- two_networks %>%
  group_by(method, hands, net1, net2) %>%
  summarize(
    n = n(),
    mean_mcc = mean(mcc),
    mean_mcc_hand = mean(mcc_hand),
    mean_mcc_hemi = mean(mcc_hemi)
  ) %>%
  ungroup() %>%
  mutate(
    net1 = ordered(net1, levels = nets),
    net2 = ordered(net2, levels = nets),
    net_x = if_else(net1 > net2, net1, net2),
    net_y = if_else(net1 > net2, net2, net1)
  )



ggplot(two_networks_summary, aes(x = net_x, y = net_y)) +
  geom_tile(aes(fill = (mean_mcc))) +
  geom_text(aes(label = round(mean_mcc, 2)), size = 2) +
  viridis::scale_fill_viridis() +
  scale_x_discrete(limits = nets) +
  scale_y_discrete(limits = nets) +
  facet_grid(cols = vars(method), rows = vars(hands)) +
  theme_bw()  +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Overall MCC")

ggplot(filter(two_networks_summary, hands == "all"), aes(x = net_x, y = net_y)) +
  geom_tile(aes(fill = mean_mcc_hand)) +
  geom_text(aes(label = round(mean_mcc_hand, 2)), size = 2) +
  viridis::scale_fill_viridis() +
  scale_x_discrete(limits = nets) +
  scale_y_discrete(limits = nets) +
  facet_grid(cols = vars(method)) +
  theme_bw()  +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Predicting handedness", fill = "MCC")

two_networks2 <- two_networks %>%
  mutate(
    LN = net1 == "LN" | net2 == "LN"
  )

ggplot(two_networks2, aes(x = size, y = mcc)) +
  geom_point(aes(color = LN)) +
  geom_smooth() +
  scale_x_log10() +
  facet_grid(cols = vars(method), rows = vars(hands)) +
  theme_bw()

two_networks3 <- two_networks2 %>%
  mutate(
    svm = method == "svm",
    lda = method == "lda",
    nn = method == "nn",
    Aud = net1 == "Aud" | net2 == "Aud",
    CO = net1 == "CO" | net2 == "CO",
    DAN = net1 == "DAN" | net2 == "DAN",
    DMN = net1 == "DMN" | net2 == "DMN",
    FPN = net1 == "FPN" | net2 == "FPN",
    LN = net1 == "LN" | net2 == "LN",
    SMN = net1 == "SMN" | net2 == "SMN",
    VIS = net1 == "VIS" | net2 == "VIS",
  )

lm1 <- lm(mcc ~ lda + svm + nn + hands + size + Aud + CO + DAN + DMN + FPN + LN +
            SMN + VIS,
          data = two_networks3)

# Get bootstrapping ====

bootstrap_files <- list.files("results/bootstrapping/",
                              pattern = ".*_bs-.*.csv",
                              full.names = TRUE)

bootstrap_results <- tibble(f = bootstrap_files) %>%
  slice_sample(prop = .1) %>%
  mutate(
    data = map(f, read_csv, show_col_types = FALSE,
               name_repair = "unique_quiet", .progress = TRUE),
    f = basename(f) %>%
      str_remove(".csv")
  )

bootstrap_big <- bootstrap_results %>%
  mutate(

    n = map_int(data, nrow),

    method = str_extract(f, "method-[^_]*"),
    # outcome = str_extract(f, "outcome-[^_]*"),
    group = str_extract(f, "test-[^_]*"),
    size = str_extract(f, "bs-[^_]*"),
    rep = str_extract(f, "rep-[^_]*"),

    # hands = replace_na(str_extract(f, "hands-[^_]*"), "all"),
    # hemis = replace_na(str_extract(f, "hemi-[^_]*"), "all"),

    across(c(method, group, size, rep), ~str_remove(.x, "^.*-")),
    across(c(size, rep), as.numeric)

  ) %>%
  unnest(data) %>%
  mutate(
    ground_hand = str_remove(ground, ".H$"),
    ground_hemi = str_remove(predicted, "^(left|right)y"),
    pred_hand = str_remove(predicted, ".H$"),
    pred_hemi = str_remove(predicted, "^(left|right)y"),
    across(c(starts_with("ground"), starts_with("pred")), as.factor),
  ) %>%
  select(-`...1`, -f) %>%
  group_by(method, group, size, rep, n)

zero_rows <- function(confusion, rows) {

  x <- confusion
  x[rows, ] <- 0

  return(x)

}

bootstrap_nested <- bootstrap_big %>%
  nest() %>%
  select(method, group, size, rep, data, n) %>%
  ungroup() %>%
  mutate(

    confusion_hemis4 = map(data,
                           ~as.data.frame.matrix(table(.x$ground, .x$predicted)),
                           .progress = TRUE),

    confusion_L = map(confusion_hemis4, ~zero_rows(.x, 3:4), .progress = TRUE),
    confusion_R = map(confusion_hemis4, ~zero_rows(.x, 1:2),.progress = TRUE),

    mcc_hemis4_overall = map_dbl(confusion_hemis4, ~mcc(confusionM = .x),
                                 .progress = TRUE),
    mcc_hemis4_L = map_dbl(confusion_L, ~mcc(confusionM = .x),
                           .progress = TRUE),
    mcc_hemis4_R = map_dbl(confusion_R, ~mcc(confusionM = .x),
                           .progress = TRUE),

  )

bootstrap_mccs <- bootstrap_nested %>%
  select(method, group, size, rep, starts_with("mcc")) %>%
  pivot_longer(starts_with("mcc"), values_to = "mcc") %>%
  mutate(
    name = str_remove(name, "mcc_hemis4_")
  )

bootstrap_cdf <- bootstrap_mccs %>%
  group_by(method, name, size) %>%
  nest() %>%
  mutate(
    ecdf = map(data, ~ecdf(.x$mcc)),
    mean = map_dbl(data, ~mean(.x$mcc)),
    sd = map_dbl(data, ~sd(.x$mcc))
  )

write_rds(bootstrap_mccs, "bootstrapped_mccs.rds")
write_rds(bootstrap_mccs, "bootstrapped_ecdfs.rds")

ggplot(bootstrap_mccs, aes(x = size, y = mcc)) +
  geom_point(alpha = 0.01, shape = 16) +
  scale_x_log10() +
  facet_grid(cols = vars(method), rows = vars(name)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Actual results ====

networks_mccs <- networks %>%
  select(outcome, method, group, hands, hemis, net1, net2, size, starts_with("mcc")) %>%
  mutate(
    within_net = TRUE
  )

ggplot(NULL) +
  geom_point(data = filter(bootstrap_mccs, name == "overall"), aes(x = size, y = mcc),
             alpha = 0.1, shape = 16) +
  geom_point(data = filter(networks_mccs, hands == "all"),
             aes(x = size, y = mcc, color = net1),
             shape = 16) +
  geom_smooth(data =  filter(bootstrap_mccs, name == "overall"), aes(x = size, y = mcc),
              alpha = 0.01, color = "black") +
  scale_x_log10() +
  facet_wrap(vars(method)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

networks_to_z <- left_join(networks_mccs, bootstrap_cdf) %>%
  mutate(
    mcc_z = map2_dbl(ecdf, mcc, ~.x(.y)),
    mcc_z2 = (mcc - mean) / sd
  )

ggplot(data = networks_to_z, aes(x = size, y = mcc_z2, color = net1)) +
  geom_point(shape = 16) +
  geom_smooth() +
  geom_hline(yintercept = 0, color = "red") +
  scale_x_log10() +
  facet_grid(cols = vars(method), rows = vars(name),
             scales = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
