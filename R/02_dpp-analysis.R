#' ---
#' title: "DPP data analysis"
#' author: Denis Agniel
#' output: github_document
#' ---

#+ r setup, warning = FALSE, message = FALSE
library(lsa)
library(here)
library(zeallot)
library(refund)
library(fda.usc)
library(tidyverse)



dpp <- read_csv(here('data', 'dpp_long_data_clean.csv'))

trt_ds <- dpp %>%
  filter(a == 1) %>%
  mutate(id = factor(id))
ctrl_ds <- dpp %>%
  filter(a == 0) %>%
  mutate(id = factor(id))

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

gam_ctrl <- mgcv::gam(x ~ s(tt) + s(tt, by = id, bs = 're'), data = ctrl_ds)
gam_trt <- mgcv::gam(x ~ s(tt) + s(tt, by = id, bs = 're'), data = trt_ds) 

smoothed_tt <- tibble(
  id = rep(unique(dpp$id), each = 101),
  tt = rep(seq(0, max(dpp$tt), length = 101), length(unique(dpp$id)))
) %>%
  inner_join(
    dpp %>% select(id, a) %>% unique
  )
smoothed_ctrl <- smoothed_tt %>%
  filter(a == 0)
smoothed_ctrl <- smoothed_ctrl %>%
  mutate(xhat = predict(gam_ctrl, newdata = smoothed_ctrl))
smoothed_trt <- smoothed_tt %>%
  filter(a == 1)
smoothed_trt <- smoothed_trt %>%
  mutate(xhat = predict(gam_trt, newdata = smoothed_trt))
smoothed_dpp <- bind_rows(
  list(smoothed_ctrl, smoothed_trt)
)
  
                  

ggplot(dpp, aes(x = tt, y = x, group = id)) +
  geom_point(alpha = 0.5) +
  geom_line(alpha = 0.1) +
  facet_wrap(~ a) +
  theme_bw() +
  ggtitle('Raw data by treatment group')

ggplot(smoothed_dpp, aes(x = tt, y = xhat, group = id)) +
  geom_line(alpha = 0.5) +
  facet_wrap(~ a) +
  theme_bw() +
  ggtitle('Smoothed data by treatment group')



obs_lin_res <-
  fit_linear_model(y_t = y_t, y_c = y_c, X_t = trt_xhat_wide, X_c = ctrl_xhat_wide) %>%
  mutate(setting = 'obs_linear')
obs_lin_res

obs_fgam_res <-
  fit_fgam(y_t = y_t, y_c = y_c, X_t = trt_xhat_wide, X_c = ctrl_xhat_wide) %>%
  mutate(setting = 'obs_fgam')
obs_fgam_res

obs_kernel_res <-
  fit_kernel_model(y_t = y_t, y_c = y_c, X_t = trt_xhat_wide, X_c = ctrl_xhat_wide) %>%
  mutate(setting = 'obs_kernel')
obs_kernel_res
