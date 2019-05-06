Bias results
================
dagniel
2019-05-03

``` r
library(knitr)
opts_chunk$set(warning = FALSE, message = FALSE, cache = FALSE, fig.width = 7, fig.height = 7)
```

``` r
library(tidyverse)
library(here)

sim_parameters <- expand.grid(
  run = 1:1000,
  n = c(50, 250, 500, 1000),
  n_i = c(5, 10, 25),
  m = c('linear', 'nonlinear'),
  s_y = 1,
  s_x = 1,
  delta = c(5, 15, 25),
  B = 0
) %>%
  filter(n == 50 | (B == 0)) %>%
  rownames_to_column('sim')

sim_res <- readRDS(here('results/linear_and_nonlinear_sims.rds'))

sim_l <- list()
for (i in 1:length(sim_res)) {
  this_sim <- sim_res[[i]]
  if (!is.null(this_sim) & class(this_sim) != 'error') {
    sim_l[[i]] <- this_sim %>%
      mutate(sim = as.character(i)) %>%
      inner_join(sim_parameters, by = 'sim')
  }
}

sim_all <- sim_l %>%
  bind_rows

unfinished_sims <- sim_parameters %>%
  anti_join(sim_all)
write_csv(unfinished_sims, 
          here('results/unfinished-20-sims.csv'))


final_sum <- sim_all %>%
  group_by(n, n_i, m, s_x, s_y, delta, estimator, type) %>%
  summarise(deltahat_s = median(delta_s, na.rm = TRUE),
            Deltahat = median(deltahat),
            Rhat = median(R, na.rm = TRUE),
            Delta = unique(delta) + 5.812196,
            delta_bias = median(delta_s - delta, na.rm = TRUE),
            deltahat_s_se = sd(delta_s, na.rm = TRUE),
            deltahat_se = sd(deltahat, na.rm = TRUE),
            Rhat_se = sd(R, na.rm = TRUE),
            R = 1-unique(delta/Delta),
            R_bias = Rhat - R,
            nsim = length(unique(run))) %>%
  ungroup %>%
  mutate(estimator = str_remove(estimator, 'delta_s_'))
write_csv(final_sum, 
          here('results/linear_and_nonlinear_bias_results.csv'))
```
