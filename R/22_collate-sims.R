#' ---
#' title: "Collect simulations from tmp"
#' output: github_document
#' ---
#' 
library(knitr)
opts_chunk$set(warning = FALSE, message = FALSE, cache = FALSE, fig.width = 7, fig.height = 7)

#'
library(dplyr)
library(here)
library(glue)
library(purrr)
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
  filter(n == 50 | (B == 0))
tmpdir <- '/n/data1/hms/dbmi/zaklab/dma12/longitudinal-surrogates/tmp/'

sim_res <- map(1:nrow(sim_parameters), function(i) {
  # browser()
  sp <- sim_parameters[i,]
  n <- sp$n
  n_i <- sp$n_i
  m <- sp$m
  s_x <- sp$s_x
  s_y <- sp$s_y
  delta <- sp$delta
  B <- sp$B
  run <- sp$run
  
  this_sim <- try(readRDS(glue('{tmpdir}/res_n{n}_ni{n_i}_m-{m}_sy{s_y}_sx{s_x}_delta{delta}_B{B}_{run}.rds')
  ))
  if (class(this_sim) == 'try-error') {
    return(NULL)
  } else return(this_sim %>%
                  mutate(n = n,
                         n_i = n_i,
                         m = m,
                         s_x = s_x,
                         s_y = s_y,
                         delta = delta,
                         B = B,
                         run = run))
})

final_sims <- sim_res %>% bind_rows
readr::write_csv(final_sims, 
                 here('results/linear_and_nonlinear_sim_results_prelim.csv'))

final_sum <- final_sims %>%
  group_by(n, n_i, m, s_x, s_y, delta, estimator, type) %>%
  summarise(deltahat_s = mean(delta_s, na.rm = TRUE),
            Deltahat = mean(deltahat),
            Rhat = mean(R, na.rm = TRUE),
            Delta = delta + 5.812196,
            R = delta/Delta,
            est_bias = mean(delta_s - delta),
            nsim = length(unique(run))) %>%
  ungroup
final_sum %>%
  sample_n(1) %>%
  select(n, n_i, m, delta) %>%
  inner_join(final_sum)
