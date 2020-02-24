devtools::install_github('denisagniel/longsurr')

library(tidyverse)
library(here)
library(zeallot)
library(longsurr)
library(refund)
library(fda.usc)
library(Rsurrogate)
library(clustermq)

select <- dplyr::select
analysis_data <- read_csv(here('data/hiv-analysis-data.csv')) %>%
  arrange(id, tt)
all_ids <- analysis_data %>%
  select(id) %>% 
  unique

c(trt_xhat_wide, ctrl_xhat_wide) %<-%
  presmooth_data(obs_data = analysis_data, 
                 options = 
                   list(plot = TRUE, 
                        # methodBwCov = 'GCV',
                        methodBwMu = 'CV',
                        methodSelectK = 'AIC',
                        useBinnedCov = FALSE,
                        verbose = TRUE,
                        usergrid = FALSE,
                        nRegGrid = 51))

smoothed_data <- as_tibble(trt_xhat_wide, rownames = 'id') %>%
  mutate(a = 1) %>%
  full_join(
    as_tibble(ctrl_xhat_wide, rownames = 'id') %>%
      mutate(a = 0)
  ) %>%
  gather(tt, x, -id, -a) %>%
  mutate(tt = as.numeric(tt))

ttl <- read_rds(here('data/time-grids-for-cv.rds'))



options(
  clustermq.defaults = list(ptn="short",
                            log_file="Rout/log%a.log",
                            time_amt = "3:00:00"
  )
)


sim_res <- Q(longsurr:::hiv_cv,
             time_list = ttl,
                  n_jobs = 500,
                  const = list(
                    all_ids = all_ids,
                    analysis_data = analysis_data,
                    smoothed_data = smoothed_data,
                    trt_xhat_wide = trt_xhat_wide
                  ),
                  pkgs=c('tidyverse',
                         'here',
                         'zeallot',
                         'longsurr',
                         'refund',
                         'fda.usc',
                         'Rsurrogate'))
write_rds(sim_res, here('results/hiv-cv-results.rds'))  
