library(tidyverse)
library(here)

bl <- read_csv(here(
  'data/data_for_denis/baseline.csv'
),
col_names = c(
  'pt_id',
  'age',
  'race_eth',
  'randomization_date',
  'site',
  'sex',
  'iv_drug_use',
  'hemophiliac',
  'homosexual',
  'screening_cd4',
  'avg_cd4',
  'base50',
  'wt_kg',
  'karnofsky_score',
  'arv_history_strat',
  'pre_arv_days',
  'zdv_30prior',
  'non_zdv_arv_prior',
  'zdv_prior175',
  'arv_history_binary',
  'symptomatic'),
col_types = cols(
  randomization_date = col_character(),
  pt_id = col_character(),
  race_eth = col_character(),
  site = col_character()
))

cd_data <- read_csv(
  here('data/data_for_denis/cd4.csv'),
  col_names = c('pt_id', 'cd4', 'date'),
  col_types = cols(
    cd4 = col_double(),
    pt_id = col_character(),
    date = col_date(format = '%d-%b-%y')
  )
)
    
    ,
  col_types = c('double', 'character', 'date')
)
