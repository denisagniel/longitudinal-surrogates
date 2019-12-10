library(dplyr)
library(clustermq)
library(here)
library(purrr)
library(readr)
library(glue)
library(ggplot2)
remotes::install_github('denisagniel/longsurr')
library(longsurr)

tmpdir <- '/n/data1/hms/dbmi/zaklab/dma12/longitudinal-surrogates/tmp/'
fs::dir_create(tmpdir)

# longsurr:::lsa_sim(n = 250, n_i = 25, m = 'nonlinear', s_y = 1, s_x = 1, delta = 15, B = 0, run = 12, tmpdir)

sim_parameters <- expand.grid(
  run = 1:1000,
  n = 250,
  n_i = c(5, 10, 25),
  m = c('linear', 'nonlinear'),
  s_y = 1,
  s_x = 1,
  delta = 25,
  B = 1000
)

options(
  clustermq.defaults = list(ptn="long",
                            log_file="Rout/log%a.log",
                            time_amt = "240:00:00"
  )
)

sim_res <- Q(longsurr:::lsa_sim, 
             n = sim_parameters$n,
             n_i = sim_parameters$n_i,
             m = sim_parameters$m,
             s_y = sim_parameters$s_y,
             s_x = sim_parameters$s_x,
             B = sim_parameters$B,
             delta = sim_parameters$delta,
             run = sim_parameters$run,
             const = list(tmpdir = tmpdir),
             n_jobs = 800,
             memory = 4000,
             fail_on_error = FALSE
)
saveRDS(sim_res, 
        here('results/boot-sims.rds'))
fs::dir_delete(tmpdir)

sessionInfo()
