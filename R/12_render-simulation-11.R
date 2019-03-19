library(here)
options(
  clustermq.defaults = list(ptn="medium",
                            log_file="Rout/log%a.log",
                            time_amt = "60:00:00"
  )
)
rmarkdown::render(here(
  'R/11_main-sim.R'
))
