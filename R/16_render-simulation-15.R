library(here)
options(
  clustermq.defaults = list(ptn="medium",
                            log_file="Rout/log%a.log",
                            time_amt = "120:00:00"
  )
)
rmarkdown::render(here(
  'R/15_large-n-sim.R'
))
