#' ---
#' title: "Collect discontinuous sims"
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
library(pbapply)
boot_params <- expand.grid(
  B = c(250),
  n = c(250, 500, 1000),
  n_i = c(5, 10, 25),
  k = 1,
  delta = c(0.1, 0.25, 0.5),
  run = 1:1000
) %>%
  filter(B == 0 | n == 250)

other_params <- expand.grid(
  B = c(0),
  n = c(250, 500, 1000),
  n_i = c(5, 10, 25),
  k = 1,
  delta = c(0.1, 0.25, 0.5),
  run = 1:1000
) %>%
  filter(B == 0 | n == 250)

tmpdir <- '/n/data1/hms/dbmi/zaklab/dma12/longitudinal-surrogates/tmp/'

boot_res <- pbapply(boot_params, 1, function(sp) {
  # browser()
  n <- sp[2]
  n_i <- sp[3]
  delta <- sp[5]
  B <- sp[1]
  run <- sp[6]
  
  this_sim <- try(readRDS(
    glue('{tmpdir}/dc-res_n{n}_ni{n_i}_delta{delta}_B{B}_{run}.rds')
  ), silent = TRUE)
  if (class(this_sim) == 'try-error') {
    # if (n_i == 25 & n == 1000) browser()
    return(NULL)
  } else return(this_sim)
})

dcb_sim <- boot_res %>% bind_rows
reg_res <- pbapply(other_params, 1, function(sp) {
  # browser()
  n <- sp[2]
  n_i <- sp[3]
  delta <- sp[5]
  B <- sp[1]
  run <- sp[6]
  
  this_sim <- try(readRDS(
    glue('{tmpdir}/dc-res_n{n}_ni{n_i}_delta{delta}_B{B}_{run}.rds')
  ), silent = TRUE)
  if (class(this_sim) == 'try-error') {
    # if (n_i == 25 & n == 1000) browser()
    return(NULL)
  } else return(this_sim)
})
dcr_sim <- reg_res %>%
  bind_rows
dc_sim <- full_join(dcb_sim, dcr_sim)

delta_0 <- dc_sim %>%
  group_by(delta) %>%
  summarise(delta_0 = mean(deltahat))



dc_sum <- dc_sim %>%
  inner_join(delta_0) %>%
  group_by(n, n_i, delta, estimator, type) %>%
  summarise(
    Delta_0 = unique(delta_0),
    Delta_s0 = unique(delta),
    R_0 = 1 - unique(Delta_s0/Delta_0),
    Deltahat = mean(deltahat, na.rm = TRUE),
    Deltahat_s = mean(delta_s, na.rm = TRUE),
    Rhat = mean(R, na.rm = TRUE),
    var_Deltahat = var(deltahat, na.rm = TRUE),
    varhat_Delta = mean(var_delta, na.rm = TRUE),
    var_Deltahat_s = var(delta_s, na.rm = TRUE),
    varhat_Delta_s = mean(var_delta_s, na.rm = TRUE),
    var_Rhat = var(R, na.rm = TRUE),
    varhat_R = mean(var_R, na.rm = TRUE),
    ci_delta_s_q_covers = 
      mean(low_q_delta_s < Delta_s0 &
             high_q_delta_s > Delta_s0, na.rm = TRUE),
    qci_halflength_delta_s = 
      mean(high_q_delta_s/2 - low_q_delta_s/2, na.rm = TRUE),
    ci_delta_s_n_covers = 
      mean(delta_s - 1.96*sqrt(var_delta_s) < Delta_s0 &
             delta_s + 1.96*sqrt(var_delta_s) > Delta_s0, na.rm = TRUE),
    ci_R_q_covers = 
      mean(low_q_R < R_0 &
             high_q_R > R_0, na.rm = TRUE),
    qci_halflength_R = 
      mean(high_q_R/2 - low_q_R/2, na.rm = TRUE),
    ci_R_n_covers = 
      mean(R - 1.96*sqrt(var_R) < R_0 &
             R + 1.96*sqrt(var_R) > R_0, na.rm = TRUE),
    nsim = n()
  ) %>%
  ungroup

dc_sum %>%
  sample_n(1) %>%
  select(n, n_i, delta, estimator) %>%
  inner_join(dc_sum) %>%
  data.frame
  

readr::write_csv(dc_sim, 
                 here('results/discontinuous-sim_all-results.csv'))
readr::write_csv(dc_sum, 
                 here('results/discontinuous-sim_summary.csv'))
