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


#' Calculate quantiles of products of two regression coefficients
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

#' Calculate quantiles of z-scores of a regression coefficient
coef_z<- function(n, fit) {
  draws <- fit$draws(n)
  effect <- colMeans(draws)
  sd <- apply(draws, 2, sd)
  effect / sd
}

#' Calculate the z-score for simulated relative abundance
p_z <- function(n, fit, trt){
  draws <- reshape_onedraw(fit$draws(n)[1,]) %>%
    rename(Taxon = "col",
           sample = "row")
  
  draw_df <- left_join(draws,
                       data.frame(
                         sample = 1:length(trt),
                         trt = trt
                       ))
  
  effect <- draw_df %>%
    group_by(Taxon, trt) %>%
    summarise(u = mean(value),
              s = sd(value),
              n = n())
  split(effect, f = effect$Taxon) %>%
    lapply(function(l) {
      s_pooled <- sum( l$s ^ 2 * (l$n - 1)) / sum(l$n - 2)
      s_pooled <- sqrt(s_pooled)
      se <- s_pooled * sqrt(sum(1 / l$n))
      (l$u[2] - l$u[1]) / se
    }) %>%
    unlist
}
