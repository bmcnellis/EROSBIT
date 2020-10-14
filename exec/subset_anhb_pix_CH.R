
library(raster)
library(snow)
source("MRLC.R")

set.seed(1)

des <- "CH"

var <- "anhb"

yrs <- c(1985:2017)

beginCluster()

subset_MRLC_pixels(desert = des, variable = var, years = yrs, type = "antecedent")

endCluster()
