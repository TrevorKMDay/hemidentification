x <- read_table(pconns_LL[1], col_names = FALSE, na = "nan")

xx <- x[rowSums(!is.na(x)) > 1, colSums(!is.na(x)) > 1]

prod(dim(x))
sum(is.na(x))
sum(is.na(x[,1]))

y <- x %>%
    mutate(
        row = row_number()
    ) %>%
    pivot_longer(-row) %>%
    mutate(
        name = as.numeric(str_remove(name, "X"))
    )

ggplot(y, aes(x = row, y = name, fill = value)) +
    geom_tile()
