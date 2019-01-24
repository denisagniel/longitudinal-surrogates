DPP data analysis
================
Denis Agniel
Thu Jan 24 12:05:41 2019

``` r
library(lsa)
library(here)
library(zeallot)
library(refund)
library(fda.usc)
library(tidyverse)



dpp <- read_csv(here('data', 'dpp_long_data_clean.csv'))

trt_ds <- dpp %>%
  filter(a == 1)
ctrl_ds <- dpp %>%
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
  presmooth_data(obs_data = dpp, n_trt = n_trt, n_ctrl = n_ctrl)
```

![](02_dpp-analysis_files/figure-markdown_github/r%20setup-1.png)![](02_dpp-analysis_files/figure-markdown_github/r%20setup-2.png)

``` r
smoothed_dpp <- as_tibble(trt_xhat_wide, rownames = 'id') %>%
  mutate(a = 1) %>%
  full_join(
    as_tibble(ctrl_xhat_wide, rownames = 'id') %>%
      mutate(a = 0)
  ) %>%
  gather(tt, x, -id, -a) %>%
  mutate(tt = as.numeric(tt))
                  

ggplot(dpp, aes(x = tt, y = x, group = id)) +
  geom_point(alpha = 0.5) +
  geom_line(alpha = 0.1) +
  facet_wrap(~ a) +
  theme_bw() +
  ggtitle('Raw data by treatment group')
```

![](02_dpp-analysis_files/figure-markdown_github/r%20setup-3.png)

``` r
ggplot(smoothed_dpp, aes(x = tt, y = x, group = id)) +
  geom_line(alpha = 0.5) +
  facet_wrap(~ a) +
  theme_bw() +
  ggtitle('Smoothed data by treatment group')
```

![](02_dpp-analysis_files/figure-markdown_github/r%20setup-4.png)

``` r
obs_lin_res <-
  fit_linear_model(y_t = y_t, y_c = y_c, X_t = trt_xhat_wide, X_c = ctrl_xhat_wide) %>%
  mutate(setting = 'obs_linear')
```

    ## Called from: fit_linear_model(y_t = y_t, y_c = y_c, X_t = trt_xhat_wide, X_c = ctrl_xhat_wide)
    ## debug at /Users/dagniel/Projects/lsa/R/fit_functional_models.R#3: n_t <- length(y_t)
    ## debug at /Users/dagniel/Projects/lsa/R/fit_functional_models.R#4: n_c <- length(y_c)
    ## debug at /Users/dagniel/Projects/lsa/R/fit_functional_models.R#5: fit_t <- pfr(y_t ~ lf(X_t, k = 30, bs = "ps"))
    ## debug at /Users/dagniel/Projects/lsa/R/fit_functional_models.R#6: fit_c <- pfr(y_c ~ lf(X_c, k = 30, bs = "ps"))
    ## debug at /Users/dagniel/Projects/lsa/R/fit_functional_models.R#12: Xbar_t <- colMeans(X_t) %>% matrix(1, ncol(X_t))
    ## debug at /Users/dagniel/Projects/lsa/R/fit_functional_models.R#13: predict_yhat_tc <- predict(fit_c, newdata = list(X_c = Xbar_t), 
    ##     type = "response", se.fit = TRUE)
    ## debug at /Users/dagniel/Projects/lsa/R/fit_functional_models.R#15: Xbar_c <- colMeans(X_c) %>% matrix(1, ncol(X_c))
    ## debug at /Users/dagniel/Projects/lsa/R/fit_functional_models.R#16: predict_yhat_ct <- predict(fit_t, newdata = list(X_t = Xbar_c), 
    ##     type = "response", se.fit = TRUE)
    ## debug at /Users/dagniel/Projects/lsa/R/fit_functional_models.R#18: data.frame(type = c("mu_t", "mu_c", "mu_st", "mu_sc", "delta", 
    ##     "delta_s", "R"), est = c(mean(y_t), mean(y_c), predict_yhat_ct$fit, 
    ##     predict_yhat_tc$fit, mean(y_t) - mean(y_c), mean(y_t) - predict_yhat_ct$fit, 
    ##     1 - (mean(y_t) - predict_yhat_ct$fit)/(mean(y_t) - mean(y_c))), 
    ##     se = c(sd(y_t)/sqrt(n_t), sd(y_c)/sqrt(n_c), predict_yhat_ct$se.fit, 
    ##         predict_yhat_tc$se.fit, sqrt(var(y_t)/n_t + var(y_c)/n_c), 
    ##         NA, NA))

``` r
obs_lin_res
```

    ##      type         est        se    setting
    ## 1    mu_t 102.7802691 0.4464607 obs_linear
    ## 2    mu_c 105.8801020 0.4434114 obs_linear
    ## 3   mu_st 102.5663359 3.7779137 obs_linear
    ## 4   mu_sc 106.3564929 3.8629109 obs_linear
    ## 5   delta  -3.0998330 0.6292383 obs_linear
    ## 6 delta_s   0.2139331        NA obs_linear
    ## 7       R   1.0690144        NA obs_linear

``` r
obs_fgam_res <-
  fit_fgam(y_t = y_t, y_c = y_c, X_t = trt_xhat_wide, X_c = ctrl_xhat_wide) %>%
  mutate(setting = 'obs_fgam')
obs_fgam_res
```

    ##    type      est         se  setting
    ## 1  mu_t 102.7803  0.4464607 obs_fgam
    ## 2  mu_c 105.8801  0.4434114 obs_fgam
    ## 3 mu_st 108.8288  3.8295274 obs_fgam
    ## 4 mu_sc 109.5746 18.3392316 obs_fgam

``` r
obs_kernel_res <-
  fit_kernel_model(y_t = y_t, y_c = y_c, X_t = trt_xhat_wide, X_c = ctrl_xhat_wide) %>%
  mutate(setting = 'obs_kernel')
obs_kernel_res
```

    ##    type      est        se    setting
    ## 1  mu_t 102.7803 0.4464607 obs_kernel
    ## 2  mu_c 105.8801 0.4434114 obs_kernel
    ## 3 mu_st      NaN       NaN obs_kernel
    ## 4 mu_sc      NaN       NaN obs_kernel
