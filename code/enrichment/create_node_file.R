library(tidyverse)

setwd("~/MyDrive/Projects/hemisphere_fingerprinting/code/enrichment")

networks0 <- read_csv(file = paste0("~/MyDrive/Projects/",
                                    "hemisphere_fingerprinting/data/",
                                    "ColeAnticevicNetPartition/",
                                    "all_parcels_to_networks.csv"),
                      show_col_types = FALSE)

colors0 <- tribble(
    ~L, ~hex,
    "Auditory", "#F98BFC",
    "CinguloOpercular", "#C350BF",
    "Default", "#FF4E31",
    "DorsalAttention", "#1EEE3A",
    "Frontoparietal", "#FDF66D",
    "Language", "#05A0A1",
    # 7f7f7f
    "OrbitoAffective", "grey50",
    "Somatomotor", "#62FAFD",
    "Visual2", "#425BFB",
  )

colors <- colors0 %>%
  mutate(
    # Convert hex to RGB, BNV expects [0, 1], not [0, 255], round for
    # readability
    rgb = map(hex, ~as_tibble(round(t(col2rgb(.x) / 255), 3))),
  ) %>%
  unnest(rgb) %>%
  select(-hex)

write_tsv(select(colors, -L), "colors.txt", col_names = FALSE)

## Nodes ====

# X, Y, Z, color, size, name

# https://neuroimaging-core-docs.readthedocs.io/en/latest/pages/atlases.html
# https://www.lead-dbs.org/helpsupport/knowledge-base/atlasesresources/cortical-atlas-parcellations-mni-space/

cogs1 <- read_table(paste0("~/Projects/glasser/mmp1.0_mniprojections/",
                            "MMP_in_MNI_corr/coords.txt"),
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

# Filter all nodes down by right joining to the the list of nodes to keep,
my_nodes <- left_join(nodes_to_keep, all_nodes,
                      by = join_by(node == label2)) %>%
  left_join(
    select(networks0, label2, L),
    by = join_by(node == label2)
  ) %>%
  group_by(L) %>%
  mutate(
    # Create individual colors for each network
    color = cur_group_id(),
    # Replace `size` column with `degree` values
    size = degree
  ) %>%
  select(X, Y, Z, color, size, node, L)

networks <- my_nodes %>%
  select(L, color) %>%
  distinct() %>%
  arrange(color)

my_nodes_final <- my_nodes %>%
  ungroup() %>%
  select(-L) %>%
  mutate(
    across(c(X, Y, Z), ~round(.x, 2))
  )

write_tsv(my_nodes_final, "enrichment_result.node", col_names = FALSE)

