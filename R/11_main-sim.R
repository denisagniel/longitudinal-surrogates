#' ---
#' title: "Main longitudinal surrogate simulation"
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
remotes::install_github('denisagniel/longsurr')


fs::dir_create(here('results'))

lsa_sim <- function(n, n_i, w, m, B, run) {

  library(dplyr)
  library(here)
  library(purrr) 
  library(zeallot)
  library(readr)
  library(glue)
  library(longsurr)
  library(refund)
  library(mgcv)
  library(Rsurrogate)
  library(tidyr)
  library(fda.usc)
  select <- dplyr::select

  print(glue('This is sim {run} for sample size {n}, number of observations {n_i}, w {w}, under model {m}, using {B} bootstrap samples.'))
    
  c(y_t, y_c, X_t, X_c, X_err_t, X_err_c, full_data, obs_data, 
    y_c_on_t, y_t_on_c, Xi_c, Xi_t) %<-%
    generate_data(n = n, n_i = n_i, Cx_0, Cx_1, w = w, sx = 0.1, 
                  sy = 0.1, grid_size = 99, nonlinear = (m == 'nl'), 
                  fpc_gen = (m == 'fpc'))
  
  select <- dplyr::select
  c(trt_xhat_wide, ctrl_xhat_wide, trt_scores, ctrl_scores) %<-%
    presmooth_data(obs_data,n = n)
  
  oracle_ests <- fit_oracle_models(y_t, y_c, X_t, X_c)
  obs_ests <- fit_obs_models(y_t, y_c, obs_data, trt_xhat_wide,
                             ctrl_xhat_wide)
  id_data <- full_data %>%
    select(id, y, a) %>%
    distinct
  if (B > 0) {
    boot_list <- map(1:B, function(b) {
      print(b)
      boot_data <- id_data %>%
        sample_frac(replace = TRUE) %>%
        mutate(old_id = id,
               id = 1:(2*n))
      by_t <- boot_data %>%
        filter(a == 1) %>%
        pull(y)
      by_c <- boot_data %>%
        filter(a == 0) %>%
        pull(y)
      boot_obs_data <- boot_data %>%
        merge(obs_data, by.x = c('old_id', 'a', 'y'), 
              by.y = c('id', 'a', 'y')) 
      b_ntrt <- length(by_t)
      b_nctrl <- length(by_c)
      
      
      c(b_trt_xhat, b_ctrl_xhat, b_trt_scores, b_ctrl_scores) %<-%
        presmooth_data(boot_obs_data,n_trt = b_ntrt, n_ctrl = b_nctrl)
      
      
      boot_obs_ests <- fit_obs_models(by_t, by_c, boot_obs_data,
                                      b_trt_xhat, b_ctrl_xhat)
      boot_obs_ests
    })
    
    boot_ests <- bind_rows(boot_list, .id = 'boot')
    # browser()
    
    boot_vars <- boot_ests %>%
      group_by(type, setting) %>%
      summarise(boot_var = var(est, na.rm = TRUE),
                low_q = quantile(est, 0.025, na.rm = TRUE),
                high_q = quantile(est, 0.975, na.rm = TRUE))
    boot_wide <- boot_ests %>%
      select(type, est, setting, boot) %>%
      spread(type, est)
    boot_sigma <- boot_wide %>%
      group_by(setting) %>%
      summarise(var_delta = var(delta, na.rm = TRUE),
                var_delta_s = var(delta_s, na.rm = TRUE),
                cov_dd = cov(delta, delta_s, use = 'p'))
    c_ds <- boot_wide %>%
      inner_join(boot_sigma) %>%
      mutate(c_alpha = (delta_s - (1-R)*delta)^2 /
               (var_delta - 2*(1-R)*cov_dd + (1-R)^2*var_delta_s))
    
    
    full_res <- full_join(oracle_ests,
                          obs_ests) %>%
      full_join(boot_vars)
  } else full_res <- full_join(oracle_ests,
                               obs_ests)
  saveRDS(full_res, here('results/scratch_res.rds'))
  full_res
  
}

# lsa_sim(n = 50, n_i = 3, m = 'nl', w = 1, B = 20)

sim_parameters <- expand.grid(
  run = 1:500,
  n = c(50, 100, 250),
  n_i = c(3, 10),
  m = c('nl'),
  w = 1,
  B = 250
)

sim_res <- Q(lsa_sim, 
             n = sim_parameters$n,
             n_i = sim_parameters$n_i,
             m = sim_parameters$m,
             w = sim_parameters$w,
             B = sim_parameters$B,
             run = sim_parameters$run,
             const = list(),
             n_jobs = 500,
             memory = 8000,
             fail_on_error = FALSE
)
saveRDS(sim_res, 
        here('results/all_sims.rds'))
sim_ds <- bind_rows(sim_res)
write_csv(sim_ds,
          here('results/all_sims.csv'))

sessionInfo()
