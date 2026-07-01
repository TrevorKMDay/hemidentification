library(tidyverse)

scores <- read_csv("results/data-hc_clf-lda_outcome-class_results.csv", 
                   show_col_types = FALSE)

folds <- read_csv("../folds_to_use2.csv")

scores2 <- scores %>%
  left_join(folds) %>%
  group_by(group) %>%
  mutate(
    across(starts_with("LD"), ~scale(.x)[, 1], .names = "{.col}_Z"),
    pred_hemi = str_extract(predicted, "[LR]H$"),
    correct = pred_hemi == hemi
  ) %>%
  arrange(desc(handedness))

average_scores <- scores2 %>%
  group_by(hemi, handedness) %>%
  summarize(
    LD1_Z = mean(LD1_Z),
    LD2_Z = mean(LD2_Z),
    LD3_Z = mean(LD3_Z)
  )

ggplot(scores2, aes(x = LD1_Z, y = LD2_Z, 
                    color = handedness, shape = handedness, fill = handedness,
                    alpha = handedness)) +
  geom_point(
    size = 2.5
  ) +
  scale_shape_manual(
    values = c("righty" = 3, "lefty" = 21)
  ) +
  scale_color_manual(
    values = c("righty" = "grey", "lefty" = "black")
  ) +
  scale_fill_manual(
    values = c("righty" = "grey", "lefty" = "pink")
  ) +
  scale_alpha_manual(
    values = c("righty" = 0.8, "lefty" = 1)
  ) +
  geom_point(
    data = average_scores, 
    shape = 22, size = 4, alpha = 1, stroke = 1.5, 
    color = "black",
    aes(fill = handedness) 
  ) +
  scale_x_continuous(breaks = c(-1, 0, 1), labels = c("-1\nLH", "0", "1\nRH")) +
  theme_bw() +
  labs(x = "LD1 (\"Hemisphere axis\")", y = "LD2", color = "Handedness",
        fill = "Handedness", shape = "Handedness", alpha = "Handedness") +
  theme(legend.position = "bottom")

ggsave("plots/LD1_LD2.png", width = 6.5, height = 4)

ggplot(scores2, aes(x = LD1_Z, y = LD3_Z, 
                    color = handedness, shape = handedness, fill = handedness,
                    alpha = handedness)) +
  geom_point(
    size = 2.5
  ) +
  scale_shape_manual(
    values = c("righty" = 3, "lefty" = 21)
  ) +
  scale_color_manual(
    values = c("righty" = "grey", "lefty" = "black")
  ) +
  scale_fill_manual(
    values = c("righty" = "grey", "lefty" = "pink")
  ) +
  scale_alpha_manual(
    values = c("righty" = 0.8, "lefty" = 1)
  ) +
  geom_point(
    data = average_scores, 
    shape = 22, size = 5, alpha = 1, stroke = 1.5, 
    color = "black",
    aes(fill = handedness) 
  ) +
  scale_x_continuous(breaks = c(-1, 0, 1), labels = c("-1\nLH", "0", "1\nRH")) +
  theme_bw() +
  labs(x = "LD1 (\"Hemisphere axis\")", y = "LD3")

ggplot(scores2, aes(x = LD2_Z, y = LD3_Z, 
                    color = handedness, shape = handedness, fill = handedness,
                    alpha = handedness)) +
  geom_point(
    size = 2.5
  ) +
  scale_shape_manual(
    values = c("righty" = 3, "lefty" = 21)
  ) +
  scale_color_manual(
    values = c("righty" = "grey", "lefty" = "black")
  ) +
  scale_fill_manual(
    values = c("righty" = "grey", "lefty" = "pink")
  ) +
  scale_alpha_manual(
    values = c("righty" = 0.8, "lefty" = 1)
  ) +
  geom_point(
    data = average_scores, 
    shape = 22, size = 4, alpha = 1, stroke = 1.5, 
    color = "black",
    aes(fill = handedness) 
  ) +
  theme_bw()

# Models

scores2_byhemi <- scores2 %>%
  group_by(hemi) %>%
  nest() %>%
  mutate(

    mean_lefty = map_dbl(data, ~mean(.x$LD1_Z[.x$handedness == "lefty"])),
    mean_righty = map_dbl(data, ~mean(.x$LD1_Z[.x$handedness == "righty"])),
    mean_diff = mean_lefty - mean_righty,

    t_test1 = map(data, ~t.test(LD1_Z ~ handedness, data = .x)),
    t_test2 = map(data, ~t.test(LD2_Z ~ handedness, data = .x))
  )
