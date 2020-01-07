ACTG 175 sim preparation
================
dagniel
2020-01-06

``` r
library(knitr)
opts_chunk$set(warning = FALSE, message = FALSE, cache = FALSE, fig.width = 7, fig.height = 7)
```

``` r
library(tidyverse)
library(here)
library(zeallot)
library(longsurr)
library(refund)
library(fda.usc)
library(Rsurrogate)
library(clustermq)
analysis_data <- read_csv(here('data/hiv-analysis-data.csv')) %>%
  arrange(id, tt)
set.seed(0)
c(trt_xhat_wide, ctrl_xhat_wide) %<-%
  presmooth_data(obs_data = analysis_data, 
                 options = 
                   list(plot = FALSE, 
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
  gather(tt, xh, -id, -a) %>%
  mutate(tt = as.numeric(tt),
         id = as.numeric(id)) %>%
  left_join(analysis_data %>%
              dplyr::select(id, tt, x)) %>%
  inner_join(analysis_data %>%
               dplyr::select(id, y) %>%
               unique)
write_rds(smoothed_data, here('results/hiv-smoothed-data.rds'))

X_1 <- trt_xhat_wide[,colnames(trt_xhat_wide) %in% colnames(ctrl_xhat_wide)]
X_0 <- ctrl_xhat_wide[,colnames(ctrl_xhat_wide) %in% colnames(trt_xhat_wide)]

fdX_t <- fdata(X_1)
fdX_c <- fdata(X_0)
k <- 3
trt_guys <- analysis_data %>%
  filter(a == 1) %>%
  dplyr::select(id, y) %>%
  unique
y_t <- trt_guys %>%
  pull(y)
control_guys <- analysis_data %>%
  filter(a == 0) %>%
  dplyr::select(id, y) %>%
  unique %>%
  mutate(X = list(X_0))
y_c <- control_guys %>%
  pull(y)
kernel_fit_1 <- fregre.np.cv(fdataobj = fdX_t, y = y_t,
                             metric = longsurr:::make_semimetric_pca(k))
kernel_fit_0 <- fregre.np.cv(fdataobj = fdX_c, y = y_c,
                             metric = longsurr:::make_semimetric_pca(k))
lin_1 <- fgam(y_t ~ lf(X_1))
lin_0 <- fgam(y_c ~ lf(X_0))
fgam_fit_1 <- fgam(y_t ~ af(X_1))
fgam_fit_0 <- fgam(y_c ~ af(X_0))
k_sigma_0 <- kernel_fit_0$sr2 %>% sqrt
k_sigma_1 <- kernel_fit_1$sr2 %>% sqrt
l_sigma_0 <- lin_0$sig2 %>% sqrt
l_sigma_1 <- lin_1$sig2 %>% sqrt
g_sigma_0 <- fgam_fit_0$sig2 %>% sqrt
g_sigma_1 <- fgam_fit_1$sig2 %>% sqrt

k_trt_guys <- trt_guys %>%
  mutate(y_1 = predict(kernel_fit_1),
         y_0 = predict(kernel_fit_0, fdX_t))
k_control_guys <- control_guys %>%
  mutate(y_0 = predict(kernel_fit_0),
         y_1 = predict(kernel_fit_1, fdX_c))

k_sim_pool <- smoothed_data %>%
  full_join(k_trt_guys %>%
              full_join(k_control_guys))

k_sim_id_data <- k_sim_pool %>%
  dplyr::select(id, a, y_1, y_0) %>%
  unique
k_sim_facts <- k_sim_id_data %>%
  summarise(Delta = mean(y_1[a == 1]) - mean(y_0[a==0]),
            Delta_S = mean(y_1[a==0]) - mean(y_0[a==0]),
            R = 1 - Delta_S/Delta)

g_trt_guys <- trt_guys %>%
  mutate(y_1 = predict(fgam_fit_1),
         y_0 = predict(fgam_fit_0, newdata = list(X_0 = X_1), type = 'response'))
g_control_guys <- control_guys %>%
  mutate(y_0 = predict(fgam_fit_0),
         y_1 = predict(fgam_fit_1, newdata = list(X_1 = X_0), type = 'response'))

g_sim_pool <- smoothed_data %>%
  full_join(g_trt_guys %>%
              full_join(g_control_guys))
g_sim_id_data <- g_sim_pool %>%
  dplyr::select(id, a, y_1, y_0) %>%
  unique
g_sim_facts <- g_sim_id_data %>%
  summarise(Delta = mean(y_1[a == 1]) - mean(y_0[a==0]),
            Delta_S = mean(y_1[a==0]) - mean(y_0[a==0]),
            R = 1 - Delta_S/Delta)

l_trt_guys <- trt_guys %>%
  mutate(y_1 = predict(lin_1),
         y_0 = predict(lin_0, newdata = list(X_0 = X_1), type = 'response'))
l_control_guys <- control_guys %>%
  mutate(y_0 = predict(lin_0),
         y_1 = predict(lin_1, newdata = list(X_1 = X_0), type = 'response'))

l_sim_pool <- smoothed_data %>%
  full_join(l_trt_guys %>%
              full_join(l_control_guys))
l_sim_id_data <- l_sim_pool %>%
  dplyr::select(id, a, y_1, y_0) %>%
  unique
l_sim_facts <- l_sim_id_data %>%
  summarise(Delta = mean(y_1[a == 1]) - mean(y_0[a==0]),
            Delta_S = mean(y_1[a==0]) - mean(y_0[a==0]),
            R = 1 - Delta_S/Delta)

sim_list <- list(k_sigma_1 = k_sigma_1,
     k_sigma_0 = k_sigma_0,
     k_sim_id_data = k_sim_id_data,
     k_sim_facts = k_sim_facts,
     g_sigma_1 = g_sigma_1,
     g_sigma_0 = g_sigma_0,
     g_sim_id_data = g_sim_id_data,
     g_sim_facts = g_sim_facts,
     l_sigma_1 = l_sigma_1,
     l_sigma_0 = l_sigma_0,
     l_sim_id_data = l_sim_id_data,
     l_sim_facts = l_sim_facts)
write_rds(sim_list, here('results/hiv-sim-tools.rds'))
```
