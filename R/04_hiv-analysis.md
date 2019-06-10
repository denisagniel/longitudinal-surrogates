ACTG 175 analysis
================
dagniel
2019-06-07

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
```

There were 460 individuals in the treatment group, and they had on average 6.4456522 longitudinal surrogate observations each.

``` r
n_ctrl <- ctrl_ds %>%
  summarise(n_ctrl = length(unique(id))) %>%
  pull(n_ctrl)
y_c <- ctrl_ds %>%
  select(id, y) %>%
  unique %>%
  pull(y)
```

There were 427 individuals in the control group, and they had on average 6.5948478 longitudinal surrogate observations each.

``` r
c(trt_xhat_wide, ctrl_xhat_wide, trt_scores, ctrl_scores) %<-%
  presmooth_data(obs_data = analysis_data, 
                 options = 
                   list(plot = TRUE, 
                        # methodBwCov = 'GCV',
                        methodBwMu = 'CV',
                        methodSelectK = 'AIC',
                        useBinnedCov = FALSE,
                        verbose = TRUE))
```

![](04_hiv-analysis_files/figure-markdown_github/unnamed-chunk-4-1.png)![](04_hiv-analysis_files/figure-markdown_github/unnamed-chunk-4-2.png)

``` r
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
```

![](04_hiv-analysis_files/figure-markdown_github/unnamed-chunk-4-3.png)

``` r
ggplot(smoothed_data, aes(x = tt, y = x, group = id)) +
  geom_line(alpha = 0.5) +
  facet_wrap(~ a) +
  theme_bw() +
  ggtitle('Smoothed data by treatment group')
```

![](04_hiv-analysis_files/figure-markdown_github/unnamed-chunk-4-4.png)

``` r
overall_treatment_effect <- mean(y_t) - mean(y_c)
```

The mean change in the treatment group was -14.976087, while the mean change in the control group was -92.0351288, for an overall treatment effect of 77.0590418.

### Results on smoothed data

``` r
delta_res <- get_delta_s(y_t = y_t,
            y_c = y_c,
            X_t = trt_xhat_wide,
            X_c = ctrl_xhat_wide) %>%
  gather(method, deltahat_s) %>%
  mutate(deltahat = overall_treatment_effect,
         R = 1 - deltahat_s/deltahat)
```

    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"

``` r
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
```

    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"

``` r
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
```

    ## # A tibble: 10 x 10
    ##    method deltahat_s deltahat     R deltahat_s_l deltahat_s_h     R_l   R_h
    ##    <chr>       <dbl>    <dbl> <dbl>        <dbl>        <dbl>   <dbl> <dbl>
    ##  1 pca2        50.4      77.1 0.345        13.8          67.5  0.173  0.816
    ##  2 pca3        46.2      77.1 0.400        10.7          60.7  0.228  0.866
    ##  3 pca4        48.3      77.1 0.373         3.49         60.6  0.212  0.955
    ##  4 pca10       44.8      77.1 0.418        -1.55         63.0  0.206  1.02 
    ##  5 fgam        21.6      77.1 0.719       -24.9         306.  -3.03   1.34 
    ##  6 kfgam       21.7      77.1 0.719       -22.0          70.3 -0.0611 1.29 
    ##  7 lin         18.9      77.1 0.755       -16.0          67.8  0.231  1.25 
    ##  8 klin         9.36     77.1 0.878       -26.0          69.8  0.124  1.33 
    ##  9 mean        56.4      77.1 0.269        38.4          68.8  0.180  0.440
    ## 10 change      59.8      77.1 0.224        17.7          85.7 -0.0274 0.749
    ## # … with 2 more variables: deltahat_s_NA <int>, R_NA <int>

### Naive approach

``` r
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
```

    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"

``` r
change_res <- R.s.estimate(sone = naive_ds %>% filter(a == 1) %>% pull(change_x),
                         szero = naive_ds %>% filter(a == 0) %>% pull(change_x),
                         yone = naive_ds %>% filter(a == 1) %>% pull(y),
                         yzero = naive_ds %>% filter(a == 0) %>% pull(y))
mean_res
```

    ## $delta
    ## [1] 77.05904
    ## 
    ## $delta.s
    ## [1] 62.88512
    ## 
    ## $R.s
    ## [1] 0.1839358

``` r
change_res
```

    ## $delta
    ## [1] 77.05904
    ## 
    ## $delta.s
    ## [1] 34.88119
    ## 
    ## $R.s
    ## [1] 0.5473447
