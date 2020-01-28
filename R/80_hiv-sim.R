devtools::install_github('denisagniel/longsurr')

library(tidyverse)
library(here)
library(zeallot)
library(longsurr)
library(refund)
library(fda.usc)
library(Rsurrogate)
library(clustermq)

analysis_data <- read_csv(here('data/hiv-analysis-data.csv'))
smoothed_data <- read_rds(here('results/hiv-smoothed-data.rds'))
sim_list <- read_rds(here('results/hiv-sim-tools.rds'))
final_list <- map(1:(length(sim_list) + 2), function(i) {
  if (i %in% 1:length(sim_list)) return(sim_list[[i]])
  else if (i == length(sim_list) + 1) smoothed_data
  else analysis_data
})
names(final_list) <- c(names(sim_list), 'smoothed_data', 'analysis_data')
sim_pars <- expand.grid(mean_fn = c('kernel', 'gam', 'linear'),
                        s = 1:1000)

options(
  clustermq.defaults = list(ptn="short",
                            log_file="Rout/log%a.log",
                            time_amt = "12:00:00"
  )
)
# attach(final_list)
# tstp <- sample_n(sim_pars, 1)
# tst <- longsurr:::hiv_sim_fn(s = tstp$s, mean_fn = tstp$mean_fn)
# tst
sim_res <- Q_rows(sim_pars, 
                  longsurr:::hiv_sim_fn,
                  n_jobs = 300,
                  export=final_list,
                  pkgs=c('tidyverse',
                         'here',
                         'zeallot',
                         'longsurr',
                         'refund',
                         'fda.usc',
                         'Rsurrogate'))
write_rds(sim_res, here('results/hiv-full-sim-results.rds'))  

