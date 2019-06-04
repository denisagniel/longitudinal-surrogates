ACTG 175 analysis
================
dagniel
2019-04-23

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
get_delta_s(y_t = y_t,
            y_c = y_c,
            X_t = trt_xhat_wide,
            X_c = ctrl_xhat_wide) %>%
  gather(method, deltahat_s) %>%
  mutate(deltahat = overall_treatment_effect,
         R = 1 - deltahat_s/deltahat)
```

    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"
    ## [1] "Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation"

    ## # A tibble: 7 x 4
    ##   method         deltahat_s deltahat     R
    ##   <chr>               <dbl>    <dbl> <dbl>
    ## 1 delta_s_k           42.0      77.1 0.456
    ## 2 delta_s_fgam        21.0      77.1 0.728
    ## 3 delta_s_kfgam        8.80     77.1 0.886
    ## 4 delta_s_lin         18.4      77.1 0.761
    ## 5 delta_s_klin        21.0      77.1 0.727
    ## 6 delta_s_mean        56.4      77.1 0.268
    ## 7 delta_s_change      60.5      77.1 0.215

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
