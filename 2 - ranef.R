# This is an example script showing how the ranef estimates were generated from the fitted models

library(unmarked)
dat <- readRDS("./output/mod_murr_54.RDS")
r <- ranef(object = dat)
saveRDS(r, "./output/mod_murr_54_ranef.RDS")
