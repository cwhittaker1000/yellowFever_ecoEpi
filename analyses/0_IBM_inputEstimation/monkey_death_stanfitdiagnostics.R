library(posterior)
source("src/R/helper.R")

diags <- check_diagnostics("data/processed/monkey_death_fit.rds")
saveRDS(diags, "data/processed/monkey_death_fit_diagnostics.rds")