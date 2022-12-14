count2prob <- function(y) {
  p <- y / rowSums(y)
  as.data.frame(p) %>%
    rowid_to_column(var = "sample") %>%
    pivot_longer(!sample, names_to = "Taxon", values_to = "value")
}

post_rel_abund <- function(p_sim) {
  p_sim%>%
    as.data.frame() %>%
    pivot_longer(everything()) %>%
    mutate(
      row = str_extract(name, "[0-9]+"),
      col = str_extract(name, ",[0-9]+"),
      row = as.integer(row),
      col = as.integer(str_remove(col, ","))
    ) %>%
    rename(sample = "row",
           Taxon = "col") %>%
    mutate(Taxon = paste0("V", Taxon)) %>%
    select(value, sample,Taxon)
}
