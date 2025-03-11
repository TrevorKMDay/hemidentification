library(tidyverse)
library(readxl)
library(here)

library(ggimage)
library(extrafont)

setwd(here("code", "hypotheses"))

# Load emoji ====

images <- tibble(f = "icons/", full.names = TRUE) %>%
  mutate(
    group = str_remove(basename(f), ".png")
  )

font_import()
loadfonts(device = "all")

# Do plotting ====

hypos <- read_xlsx(here("data", "hypotheses.xlsx")) %>%
  filter(
    !is.na(hypo)
  )%>%
  rowwise() %>%
  mutate(
    same_handedness = case_when(
      sinistrality == "R" & str_detect(target, "^d") ~ "same",
      sinistrality == "L" & str_detect(target, "^s") ~ "same",
      TRUE ~ "diff"
    ),
    same_hemisphere = case_when(
      hemi == "RH" & str_detect(target, "RH$") ~ "same",
      hemi == "LH" & str_detect(target, "LH$") ~ "same",
      TRUE ~ "diff"
    ),
    group = paste(same_hemisphere, same_handedness, sep = "_"),
    score = score * runif(1, min = 0.95, max = 1.05),
    sinistrality = paste("Hand:", sinistrality)
  ) %>%
  left_join(images, by = join_by(group))

# ggplot(hypos, aes(y = hemi, x = score, shape = same_handedness,
#                   fill = same_hemisphere)) +
#   geom_point(size = 7.5, alpha = 0.75, stroke = 1.25) +
#   scale_x_continuous(limits = c(0, 1), breaks = NULL) +
#   scale_shape_manual(values = c(23, 21)) +
#   scale_fill_manual(values = c("darkred", "springgreen4")) +
#   theme_minimal() +
#   theme(panel.grid.minor = element_blank(),
#         strip.placement = "outside",
#         legend.position = "none",
#         text = element_text(size = 25)) +
#   labs(y = NULL, x = "Accuracy", shape = "Handedness", fill = "Hemisphere") +
#   facet_grid(rows = vars(sinistrality), switch = "y")

my_theme <- list(
    scale_x_continuous(limits = c(0, 1), breaks = NULL),
    theme_minimal(),
    theme(panel.grid.minor = element_blank(),
          strip.placement = "outside",
          legend.position = "none",
          text = element_text(size = 15, family = "Helvetica")),
    labs(y = NULL, x = "f(Assignment)"),
    facet_grid(rows = vars(sinistrality), switch = "y")
  )

H1 <- filter(hypos, hypo == "H1")

png("h1.png", width = 5, height = 4, units = "in", res = 300)

ggplot(H1, aes(y = hemi, x = score,
                  shape = interaction(same_handedness, same_hemisphere))) +
  geom_image(aes(image = f), size = 0.2) +
  my_theme +
  labs(title = "H1: Good performance")

dev.off()

H2 <- filter(hypos, hypo == "H2")

png("h2.png", width = 5, height = 4, units = "in", res = 300)

ggplot(H2, aes(y = hemi, x = score,
               shape = interaction(same_handedness, same_hemisphere))) +
  geom_image(aes(image = f), size = 0.2) +
  my_theme+
  labs(title = "H2: Good chirality performance")

dev.off()

H3 <- filter(hypos, hypo == "H3")

png("h3.png", width = 5, height = 4, units = "in", res = 300)

ggplot(H3, aes(y = hemi, x = score,
               shape = interaction(same_handedness, same_hemisphere))) +
  geom_image(aes(image = f), size = 0.2) +
  my_theme +
  labs(title = "H3: Good handedness performance")

dev.off()

