#' ---
#' title: "Collect simulations from tmp"
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
  run = 1:500,
  n = c(50, 100, 250),
  n_i = c(3, 10),
  m = c('nl'),
  w = 1,
  B = 250
)

all_sims <- readRDS(here('results/all_sims.rds'))

sim_res <- map(1:nrow(sim_parameters), function(i) {
  
  sp <- sim_parameters[i,]
  n <- sp$n
  n_i <- sp$n_i
  m <- sp$m
  w <- sp$w
  B <- sp$B
  run <- sp$run
  
  this_sim <- all_sims[[i]]
  if (class(this_sim) == 'error') {
    return(NULL)
  } else return(this_sim %>%
                  mutate(n = n,
                         n_i = n_i,
                         m = m,
                         w = w,
                         B = B,
                         run = run))
})

final_sims <- sim_res %>% bind_rows
readr::write_csv(final_sims, 
          here('results/all_simulation_results.csv'))
fs::dir_delete(here('results/tmp'))
