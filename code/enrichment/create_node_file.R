library(tidyverse)

setwd("~/Projects/glasser/mmp1.0_mniprojections/")

# andrew <- read_tsv("~/Projects/glasser/mmp1.0_mniprojections/example_node_file.node",
#                    col_names = FALSE, show_col_types = FALSE)
#
# range(andrew$X1)
# range(andrew$X2)
# range(andrew$X3)

networks0 <- read_csv(file = paste0("~/MyDrive/Projects/",
                                    "hemisphere_fingerprinting/data/",
                                    "ColeAnticevicNetPartition/",
                                    "all_parcels_to_networks.csv"),
                      show_col_types = FALSE)

## Nodes ====

# X, Y, Z, color, size, name

# https://neuroimaging-core-docs.readthedocs.io/en/latest/pages/atlases.html
# https://www.lead-dbs.org/helpsupport/knowledge-base/atlasesresources/cortical-atlas-parcellations-mni-space/

cogs1 <- read_table("~/Projects/glasser/mmp1.0_mniprojections/MMP_in_MNI_corr/coords.txt",
                   col_names = c("index", "X", "Y", "Z"),
                   show_col_types = FALSE) %>%
  bind_cols(networks0) %>%
  select(-c(index, L, R))

range(cogs1$X)
range(cogs1$Y)
range(cogs1$Z)

# cogs2 <- read_table("~/Projects/glasser/mmp1.0_mniprojections/HCPMMP1_on_MNI152_ICBM2009a_nlin_hd//coords.txt",
#                     col_names = c("index", "X", "Y", "Z"),
#                     show_col_types = FALSE) %>%
#   bind_cols(networks0) %>%
#   select(-c(index, L, R))
#
# range(cogs2$X)
# range(cogs2$Y)
# range(cogs2$Z)

all_nodes <- cogs1 %>%
  mutate(
    color = 1,
    size = 1,
  ) %>%
  select(X, Y, Z, color, size, label2)

write_tsv(all_nodes, "glasser_all.node", col_names = FALSE)

nodes_to_keep <- read_tsv("enrichment_nodes_to_keep.tsv",
                          show_col_types = FALSE)

my_nodes <- all_nodes %>%
  filter(
    label2 %in% nodes_to_keep$node
  ) %>%
  left_join(
    select(networks0, label2, L)
  ) %>%
  group_by(L) %>%
  mutate(
    color = cur_group_id()
  )

networks <- my_nodes %>%
  select(L, color) %>%
  distinct() %>%
  arrange(color)

my_nodes_final <- my_nodes %>%
  ungroup() %>%
  select(-L)

write_tsv(my_nodes_final, "enrichment_result.node", col_names = FALSE)

