#' ---
#' title: "ACTG 175 analysis"
#' output: github_document
#' ---
#' 
library(knitr)
opts_chunk$set(warning = FALSE, message = FALSE, cache = FALSE, fig.width = 7, fig.height = 7)

#'
library(tidyverse)
library(here)
library(zeallot)
library(longsurr)
library(refund)
library(fda.usc)
library(Rsurrogate)
select <- dplyr::select
analysis_data <- read_csv(
          here('data/hiv-analysis-data.csv'))

trt_ds <- analysis_data  %>%
  filter(a == 1)
ctrl_ds <- analysis_data %>%
  filter(a == 0)

n_trt <- trt_ds %>%
  summarise(n_trt = length(unique(id))) %>%
  pull(n_trt)
y_t <- trt_ds %>%
  select(id, y) %>%
  unique %>%
  pull(y)

#' There were `r n_trt` individuals in the treatment group, and they had on average `r nrow(trt_ds)/n_trt` longitudinal surrogate observations each. 

n_ctrl <- ctrl_ds %>%
  summarise(n_ctrl = length(unique(id))) %>%
  pull(n_ctrl)
y_c <- ctrl_ds %>%
  select(id, y) %>%
  unique %>%
  pull(y)
#' There were `r n_ctrl` individuals in the control group, and they had on average `r nrow(ctrl_ds)/n_ctrl` longitudinal surrogate observations each. 
c(trt_xhat_wide, ctrl_xhat_wide, trt_scores, ctrl_scores) %<-%
  presmooth_data(obs_data = analysis_data, 
                 options = 
                   list(plot = TRUE, 
                        # methodBwCov = 'GCV',
                        methodBwMu = 'CV',
                        methodSelectK = 'AIC',
                        useBinnedCov = FALSE,
                        verbose = TRUE))

smoothed_data <- as_tibble(trt_xhat_wide, rownames = 'id') %>%
  mutate(a = 1) %>%
  full_join(
    as_tibble(ctrl_xhat_wide, rownames = 'id') %>%
      mutate(a = 0)
  ) %>%
  gather(tt, x, -id, -a) %>%
  mutate(tt = as.numeric(tt))


ggplot(analysis_data, aes(x = tt, y = x, group = id)) +
  geom_point(alpha = 0.5) +
  geom_line(alpha = 0.1) +
  facet_wrap(~ a) +
  theme_bw() +
  ggtitle('Raw data by treatment group')

ggplot(smoothed_data, aes(x = tt, y = x, group = id)) +
  geom_line(alpha = 0.5) +
  facet_wrap(~ a) +
  theme_bw() +
  ggtitle('Smoothed data by treatment group')

overall_treatment_effect <- mean(y_t) - mean(y_c)

#' The mean change in the treatment group was `r mean(y_t)`, while the mean change in the control group was `r mean(y_c)`, for an overall treatment effect of `r overall_treatment_effect`. 
#'
#'
#' ### Results on smoothed data
delta_res <- get_delta_s(y_t = y_t,
            y_c = y_c,
            X_t = trt_xhat_wide,
            X_c = ctrl_xhat_wide) %>%
  gather(method, deltahat_s) %>%
  mutate(deltahat = overall_treatment_effect,
         R = 1 - deltahat_s/deltahat)

id_data <- analysis_data %>%
  select(id, a, y) %>%
  distinct

boot_l <- map(1:200, function(b) {
  boot_data <- id_data %>%
    sample_frac(replace = TRUE)
  boot_data <- boot_data %>%
    arrange(id) %>%
    mutate(old_id = id,
           id = 1:nrow(boot_data))
  boot_obs_data <- boot_data %>%
    merge(analysis_data, by.x = c('old_id', 'a', 'y'), 
          by.y = c('id', 'a', 'y')) %>%
    arrange(id, tt)
  
  trt_ds <- boot_obs_data  %>%
    filter(a == 1)
  ctrl_ds <- boot_obs_data %>%
    filter(a == 0)
  
  n_trt <- trt_ds %>%
    summarise(n_trt = length(unique(id))) %>%
    pull(n_trt)
  y_t <- trt_ds %>%
    select(id, y) %>%
    unique %>%
    pull(y)
  
  n_ctrl <- ctrl_ds %>%
    summarise(n_ctrl = length(unique(id))) %>%
    pull(n_ctrl)
  y_c <- ctrl_ds %>%
    select(id, y) %>%
    unique %>%
    pull(y)

  c(trt_xhat_wide, ctrl_xhat_wide, trt_scores, ctrl_scores) %<-%
    presmooth_data(obs_data = boot_obs_data, 
                   options = 
                     list(plot = FALSE, 
                          # methodBwCov = 'GCV',
                          methodBwMu = 'CV',
                          methodSelectK = 'AIC',
                          useBinnedCov = FALSE,
                          verbose = TRUE))
  
  overall_treatment_effect <- mean(y_t) - mean(y_c)
  get_delta_s(y_t = y_t,
              y_c = y_c,
              X_t = trt_xhat_wide,
              X_c = ctrl_xhat_wide) %>%
    gather(method, deltahat_s) %>%
    mutate(deltahat = overall_treatment_effect,
           R = 1 - deltahat_s/deltahat,
           boot = b)
})
boot_res <- boot_l %>%
  bind_rows %>%
  group_by(method) %>%
  summarise(deltahat_s_l = quantile(deltahat_s, 0.025, na.rm = TRUE),
            deltahat_s_h = quantile(deltahat_s, 0.975, na.rm = TRUE),
            R_l = quantile(R, 0.025, na.rm = TRUE),
            R_h = quantile(R, 0.975, na.rm = TRUE),
            deltahat_s_NA = sum(is.na(deltahat_s)),
            R_NA = sum(is.na(R)))
delta_res %>%
  inner_join(boot_res)

#' ### Naive approach

naive_ds <- analysis_data %>%
  group_by(id, y, a) %>%
  summarise(
    mean_x = mean(x),
    change_x = x[tt == max(tt)] - x[tt == min(tt)]
  )
  
mean_res <- R.s.estimate(sone = naive_ds %>% filter(a == 1) %>% pull(mean_x),
                         szero = naive_ds %>% filter(a == 0) %>% pull(mean_x),
                         yone = naive_ds %>% filter(a == 1) %>% pull(y),
                         yzero = naive_ds %>% filter(a == 0) %>% pull(y))
change_res <- R.s.estimate(sone = naive_ds %>% filter(a == 1) %>% pull(change_x),
                         szero = naive_ds %>% filter(a == 0) %>% pull(change_x),
                         yone = naive_ds %>% filter(a == 1) %>% pull(y),
                         yzero = naive_ds %>% filter(a == 0) %>% pull(y))
mean_res
change_res
