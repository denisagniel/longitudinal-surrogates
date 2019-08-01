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


# tmp <- longsurr::dc_sim(n = 500, n_i = 5, k = 1, s_y = 0.5, s_x = 0.5, delta = 0.5, B = 0, run = -1, tmpdir)
# tmp

sim_parameters <- expand.grid(
  run = 1:3,
  n = c(250, 500, 1000),
  n_i = c(5, 10, 25),
  k = 1,
  delta = c(0.1, 0.25, 0.5),
  B = c(0, 250)
) %>%
  filter(B == 0 | n == 250)

options(
  clustermq.defaults = list(ptn="short",
                            log_file="Rout/log%a.log",
                            time_amt = "12:00:00"
  )
)

# sim_res <- Q(dc_sim, 
#              n = sim_parameters$n,
#              n_i = sim_parameters$n_i,
#              k = sim_parameters$k,
#              B = sim_parameters$B,
#              delta = sim_parameters$delta,
#              run = sim_parameters$run,
#              const = list(tmpdir = tmpdir),
#              n_jobs = 500,
#              # memory = 2000,
#              fail_on_error = FALSE
# )
sim_res <- Q_rows(sim_parameters, fun = dc_sim,
                  const = list(tmpdir = tmpdir),
                  n_jobs = 50,
                  # memory = 2000,
                  fail_on_error = FALSE)
saveRDS(sim_res, 
        here('results/discontinuous_sims.rds'))
fs::dir_delete(tmpdir)

sessionInfo()
