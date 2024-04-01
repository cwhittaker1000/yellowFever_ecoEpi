library(posterior)
source("src/R/helper.R")

diags <- check_diagnostics("data/processed/observation_time_fit.rds")
saveRDS(diags, "data/processed/observation_time_fit_diagnostics.rds")