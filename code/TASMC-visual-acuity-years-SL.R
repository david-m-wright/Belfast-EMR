# Superlearner models to predict visual acuity at end of each treatment year

source(rprojroot::find_rstudio_root_file("code/TASMC-injection-demand.R"))
set.seed(1983)
# Individual models for each year

# Create sl3 tasks for each year of injections
va_task_yr1_6 <- make_sl3_Task(
  data = filter(fluid_6months_imp_demand_va_yr1, years_observed >=1),
  covariates =  noa_covariates_0_6,
  outcome = "va_logmar_yr_1",
  outcome_type = "continuous",
  id = "eye_id"
)

va_task_yr2_6 <- make_sl3_Task(
  data = filter(fluid_6months_imp_demand_va_yr2, years_observed >=2),
  covariates =  noa_covariates_0_6,
  outcome = "va_logmar_yr_2",
  outcome_type = "continuous",
  id = "eye_id"
)

va_task_yr3_6 <- make_sl3_Task(
  data = filter(fluid_6months_imp_demand_va_yr3, years_observed >=3),
  covariates =  noa_covariates_0_6,
  outcome = "va_logmar_yr_3",
  outcome_type = "continuous",
  id = "eye_id"
)

# List the demand tasks
va_tasks_yrs <- list(va_task_yr1_6, va_task_yr2_6, va_task_yr3_6) 


# Load the fitted model
sl_fit_va_yrs_6 <- read_rds(find_rstudio_root_file("SL-models", "sl_fit_va_yrs_6.RDS"))

# Extract predictions
sl_preds_va_yrs_6 <- pmap(list(va_tasks_yrs,
                               sl_fit_va_yrs_6,
                               fluid_6months_imp_demand_va_yrs),
                          ~extract_sl_continuous(..1, ..2, ..3, eye_demand, round = FALSE)) %>%
  set_names(c("va_yr_1", "va_yr_2", "va_yr_3")) %>%
  bind_rows(.id = "outcome")


