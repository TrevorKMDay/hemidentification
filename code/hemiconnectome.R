library(tidyverse)

setwd("o")

pconn <- read_delim("data/test.pconn.txt", delim = "\t", show_col_types = FALSE,
                    col_names = FALSE) 

parcels <- read_table("data/parcels.txt", skip = 1, show_col_types = FALSE,
                      col_names = c("parcel", "index", "label")) %>%
  mutate(
    index = as.numeric(str_remove(index, ":")),
    label = str_remove(label, "_ROI")
  )

pconn2 <- pconn 
colnames(pconn2) <- parcels$label

# Remove subcortical
pconn3 <- pconn2[1:360, 1:360]
pconn3[upper.tri(pconn3, diag = TRUE)] <- NA

pconn4 <- pconn3 %>%
  mutate(
    var1 = parcels$label[1:360]
  ) %>%
  pivot_longer(-var1, names_to = "var2", values_to = "z") %>%
  na.omit() %>%
  rowwise() %>%
  mutate(
    r = round(psych::fisherz2r(z), 3),
    index1 = which(parcels$label == var1),
    index2 = which(parcels$label == var2),
    
    hemi = case_when(
      str_detect(var1, "^R_") & str_detect(var2, "^R_") ~ "R_R",
      str_detect(var1, "^L_") & str_detect(var2, "^L_") ~ "L_L",
      TRUE ~ "callosal"
    )
    
  ) 

labels <- c(seq(0, 180, 60), seq(60, 180, 60))

png("hemiconnectome.png", width = 4.5, height = 3.5, units = "in", res = 300)

ggplot(pconn4, aes(x = index1, y = index2, fill = r)) +
  geom_tile() +
  scale_x_reverse(breaks = seq(0, 360, by = 60), labels = labels) +
  scale_y_continuous(breaks = seq(0, 360, by = 60), labels = labels) +
  scale_fill_gradient2(limits = c(-1, 1)) +
  coord_equal() +
  theme_bw()  +
  labs(x = "ROI 1", y = "ROI 2") +
  annotate("rect", xmin = 360, xmax = 180, ymin = 180, ymax = 360, 
            color = "red", fill = NA, size = 1) +
  annotate("rect", xmin = 0, xmax = 180, ymin = 180, ymax = 0, 
            color = "blue", fill = NA, size = 1) +
  annotate("rect", xmin = 360, xmax = 180, ymin = 180, ymax = 0, 
            color = "purple", fill = NA, size = 1)

dev.off()

lm(r ~ hemi, data = pconn4) %>%
  summary()

ggplot(pconn4, aes(x = hemi, y = r)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(NA, 1)) +
  theme_bw()
