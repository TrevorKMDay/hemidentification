library(tidyverse)
library(readxl)
library(scales)

setwd("H:/My Drive/Projects/hemisphere_fingerprinting/code")

fconn <- read_xlsx("fake_hemiconn.xlsx")

fconn2 <- fconn %>%
  pivot_longer(-ROI1, names_to = "ROI2", values_to = "r") %>%
  na.omit() %>%
  mutate(
    r = if_else(str_detect(ROI1, "LH") & str_detect(ROI2, "RH"), r - 0.5, r),
    r_noisy = jitter(r, amount = 0.2) %>%
      replace(., . > 1, 1),
    across(starts_with("ROI"), ~str_replace(.x, "_", " "))
  )

png("fake_conn_8x8.png", width = 4.5, height = 3.5, units = "in", res = 300)

ggplot(fconn2, aes(y = ROI1, x = ROI2, fill = r_noisy)) +
  geom_tile() +
  scale_y_discrete(limits = rev) +
  scale_fill_gradient2(limits = c(-1, 1), low = muted("blue"),
                       high = muted("red")) +
  coord_equal() +
  theme_bw() +
  labs(x = "ROI 1", y = "ROI 2", fill = "r")

dev.off()
