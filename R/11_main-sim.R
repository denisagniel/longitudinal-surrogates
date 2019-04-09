#' ---
#' title: "Main longitudinal surrogate simulation, nonlinear"
#' output: github_document
#' ---
#' 
library(knitr)
opts_chunk$set(warning = FALSE, message = FALSE, cache = FALSE, fig.width = 7, fig.height = 7)

#'
library(dplyr)
library(clustermq)
library(here)
library(purrr)
library(readr)
library(glue)
library(ggplot2)
library(longsurr)
remotes::install_github('denisagniel/longsurr')


fs::dir_create(here('results'))
fs::dir_create(here('results/tmp'))


longsurr::lsa_sim(n = 250, n_i = 3, k = 1, s_y = 0.5, s_x = 0.5, delta = 0.25, B = 3, run = 0)

sim_parameters <- expand.grid(
  run = 1:500,
  n = c(50, 100, 250, 500, 1000),
  n_i = c(3, 5, 10),
  k = 1,
  s_y = 0.5,
  s_x = 0.5,
  delta = c(0.1, 0.25, 0.5),
  B = c(0, 250)
) %>%
  filter(n == 50 | (B == 0))

sim_res <- Q(lsa_sim, 
             n = sim_parameters$n,
             n_i = sim_parameters$n_i,
             k = sim_parameters$k,
             s_y = sim_parameters$s_y,
             s_x = sim_parameters$s_x,
             B = sim_parameters$B,
             delta = sim_parameters$delta,
             run = sim_parameters$run,
             const = list(),
             n_jobs = 500,
             memory = 1000,
             fail_on_error = FALSE
)
saveRDS(sim_res, 
        here('results/nonlinear_sims.rds'))
sim_ds <- bind_rows(sim_res)
write_csv(sim_ds,
          here('results/nonlinear_sims.csv'))
fs::dir_delete(here('results/tmp'))

sessionInfo()
