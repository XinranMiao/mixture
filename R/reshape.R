count2prob <- function(y) {
  p <- y / rowSums(y)
  as.data.frame(p) %>%
    rowid_to_column(var = "sample") %>%
    pivot_longer(!sample, names_to = "Taxon", values_to = "value")
}

post_rel_abund <- function(p_sim) {
  p_sim %>%
    reshape_onedraw() %>%
    rename(sample = "row",
           Taxon = "col") %>%
    mutate(Taxon = paste0("V", Taxon)) %>%
    select(value, sample,Taxon)
}


reshape_onedraw <- function(x) {
  x %>%
    as.data.frame() %>%
    pivot_longer(everything()) %>%
    mutate(
      row = str_extract(name, "[0-9]+"),
      col = str_extract(name, ",[0-9]+"),
      row = as.integer(row),
      col = as.integer(str_remove(col, ","))
    )
}



#' Extract quantiles of regression coefficients
coef_quantile <- function(n, fit) {
  fit$draws(n) %>%
    apply(2, quantile, c(0.05, 0.5, 0.95)
    ) %>%
    t() %>%
    as.data.frame() %>%
    mutate(name = n) %>%
    rename(lwr = `5%`,
           upr = `95%`,
           estimate = `50%`) %>%
    rowid_to_column(var = "Taxon") %>%
    mutate(Taxon = str_extract(Taxon, "[0-9]+"),
           Taxon = str_c("V", Taxon))
}

coef_prod_quantile <- function(n1, n2, fit) {
  (fit$draws(n1) * fit$draws(n2)) %>%
    apply(2, quantile, c(0.05, 0.5, 0.95)
    ) %>%
    t() %>%
    as.data.frame() %>%
    mutate(name = str_c(n1, ".", n2)) %>%
    rename(lwr = `5%`,
           upr = `95%`,
           estimate = `50%`) %>%
    rowid_to_column(var = "Taxon") %>%
    mutate(Taxon = str_extract(Taxon, "[0-9]+"),
           Taxon = str_c("V", Taxon))
}
