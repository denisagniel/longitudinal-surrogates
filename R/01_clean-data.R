library(tidyverse)
library(lsa)
library(here)
library(janitor)

dpp_data <- read_csv(here('data', 'dpp.data.csv')) %>%
  clean_names

### make exclusions to for now focus on placebo vs lifestyle

dpp_lp <- dpp_data %>%
  filter(treatment_group %in% c('Lifestyle', 'Placebo')) %>%
  mutate(a = if_else(treatment_group == 'Lifestyle', 1, 0))

#-----
## reshape the data
#--------

dpp_lp_long_t <-
  dpp_lp %>%
  select(release_id, matches('HBA1_._time')) %>%
  gather(time_n, tt, -release_id) %>%
  separate(time_n, into = c('measure', 'time_n', 'tv')) %>%
  select(-tv)

dpp_lp_long_v <-
  dpp_lp %>%
  select(release_id, matches('HBA1_._value')) %>%
  gather(val_n, x, -release_id) %>%
  separate(val_n, into = c('measure', 'time_n', 'tv')) %>%
  select(-tv)

dpp_lp_long <- dpp_lp_long_t %>%
  inner_join(dpp_lp_long_v) %>% 
  inner_join(
    dpp_lp %>%
      transmute(release_id, 
                id = as.numeric(as.factor(release_id)),
                a,
             y = glucose_2yr_value,
             y_diff = glucose_2yr_value - glucose_baseline_value)
) %>%
  filter(!is.na(tt))

write_csv(dpp_lp_long, here('data', 'dpp_long_data_clean.csv'))
