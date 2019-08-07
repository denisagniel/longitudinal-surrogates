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
sim_parameters <- expand.grid(
  run = 1:1000,
  n = c(250, 500, 1000),
  n_i = c(5, 10, 25),
  k = 1,
  delta = c(0.1, 0.25, 0.5),
  B = 0
)

dc_sim <- readRDS(here("results/discontinuous_sims.rds")) %>%
  bind_rows 

delta_0 <- dc_sim %>%
  group_by(delta) %>%
  summarise(delta_0 = mean(deltahat))
  
dc_res <- dc_sim %>%
  group_by(n, n_i, delta, type, estimator) %>%
  summarise(deltahat_s = mean(delta_s, na.rm = TRUE),
            deltahat = mean(deltahat, na.rm = TRUE),
            Rhat = mean(R, na.rm = TRUE),
            var_deltahat_s = var(delta_s, na.rm = TRUE),
            var_deltahat = var(deltahat,na.rm = TRUE),
            var_R = var(R, na.rm = TRUE)) %>%
  inner_join(delta_0) %>%
  mutate(R_0 = 1 - delta/delta_0,
         delta_s_bias = deltahat_s - delta,
         delta_bias = deltahat - delta_0,
         R_bias = Rhat - R_0)

dc_res %>% group_by(estimator, type) %>% 
  summarise_all(~mean(.)) %>% arrange(R_bias)
readr::write_csv(final_sims, 
                 here('results/good-sim-results-cleaned-up.csv'))
