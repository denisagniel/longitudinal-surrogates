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
# 
# longsurr:::hiv_cv(s = 0, 
#                   time_list = sample(ttl, 3),
#                   all_ids = all_ids,
#                   analysis_data = analysis_data,
#                   smoothed_data = smoothed_data,
#                   trt_xhat_wide = trt_xhat_wide
#                   )


options(
  clustermq.defaults = list(ptn="short",
                            log_file="Rout/process-%a.log",
                            time_amt = "3:00:00"
  )
)
# 
# sim_res <- Q(longsurr:::hiv_cv,
#              s = 2,
#              n_jobs = 2,
#              const = list(
#                time_list = sample(ttl, 3),
#                all_ids = all_ids,
#                analysis_data = analysis_data,
#                smoothed_data = smoothed_data,
#                trt_xhat_wide = trt_xhat_wide
#              ),
#              pkgs=c('tidyverse',
#                     'here',
#                     'zeallot',
#                     'longsurr',
#                     'refund',
#                     'fda.usc',
#                     'Rsurrogate'))


sim_res <- Q(longsurr:::hiv_cv,
             s = 1:500,
                  n_jobs = 500,
                  const = list(
                    time_list = ttl,
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
                         'Rsurrogate'),
             fail_on_error = FALSE)
write_rds(sim_res, here('results/hiv-full-cv-results.rds'))  
