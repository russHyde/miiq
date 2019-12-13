message("Setting up Development env for {miiq}")

Sys.setenv(R_REMOTES_UPGRADE = "never")
library(devtools)
library(testthat)
library(goodpractice)
library(lintr)
library(styler)

gp_covr <- function(...) gp(..., checks = "covr")

