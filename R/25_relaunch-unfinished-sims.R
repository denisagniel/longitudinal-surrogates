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


# longsurr::lsa_sim(n = 50, n_i = 10, m = 'nonlinear', s_y = 1, s_x = 1, delta = 15, B = 0, run = 12, tmpdir)

sim_parameters <- read_csv(here('results/unfinished-20-sims.csv'))

options(
  clustermq.defaults = list(ptn="short",
                            log_file="Rout/log%a.log",
                            time_amt = "12:00:00"
  )
)

sim_res <- Q(lsa_sim, 
             n = sim_parameters$n,
             n_i = sim_parameters$n_i,
             m = sim_parameters$m,
             s_y = sim_parameters$s_y,
             s_x = sim_parameters$s_x,
             B = sim_parameters$B,
             delta = sim_parameters$delta,
             run = sim_parameters$run,
             const = list(tmpdir = tmpdir),
             n_jobs = 100,
             memory = 2000,
             fail_on_error = FALSE
)
saveRDS(sim_res, 
        here('results/relaunched_linear-and-nonlinear-sims.rds'))
fs::dir_delete(tmpdir)

sessionInfo()
