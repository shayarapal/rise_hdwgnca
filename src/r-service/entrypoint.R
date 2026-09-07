library(plumber)
library(future)
library(promises)

# Two workers: one for /test-soft-powers, one for /analyze running simultaneously
plan(multisession, workers = 2)

port <- as.integer(Sys.getenv("PORT", "8100"))

pr("plumber.R") |>
  pr_run(host = "0.0.0.0", port = port)
