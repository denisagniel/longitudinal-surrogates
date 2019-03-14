library(here)
options(
  clustermq.defaults = list(ptn="medium",
                            log_file="Rout/log%a.log",
                            time_amt = "36:00:00"
  )
)
rmarkdown::render(here(
  'R/11_main-sim.R'
))
