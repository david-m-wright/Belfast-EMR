# Superlearner models to predict visual acuity at end of each treatment year

source(rprojroot::find_rstudio_root_file("code/BIRAX-injection-demand.R"))
set.seed(1983)

# Define tasks
# These operate on different sized source datasets as some eyes had no VA measurements in years 2 and 3
va_task_yr1_6 <- make_sl3_Task(
  data = fluid_6months_imp_demand_va_yr1,
  covariates = noa_covariates_0_6,
  outcome = "va_logmar_yr_1",
  outcome_type = "continuous",
  id = "PatientID",
  folds = n_folds
)

va_task_yr2_6 <- make_sl3_Task(
  data = fluid_6months_imp_demand_va_yr2,
  covariates = noa_covariates_0_6,
  outcome = "va_logmar_yr_2",
  outcome_type = "continuous",
  id = "PatientID",
  folds = n_folds
)

va_task_yr3_6 <- make_sl3_Task(
  data = fluid_6months_imp_demand_va_yr3,
  covariates = noa_covariates_0_6,
  outcome = "va_logmar_yr_3",
  outcome_type = "continuous",
  id = "PatientID",
  folds = n_folds
)

# List the tasks
va_tasks_yrs <- list(va_task_yr1_6, va_task_yr2_6, va_task_yr3_6)

# Fit the superlearner
sl_fit_va_yrs_6 <- map(va_tasks_yrs, ~sl$train(.))

# Calculate the cross validated risk
sl_fit_va_yrs_6_cv_risk <- map(sl_fit_va_yrs_6, ~.$cv_risk(eval_fun = loss_squared_error))

# Extract predictions
sl_preds_va_yrs_6 <- pmap(list(va_tasks_yrs,
                               sl_fit_va_yrs_6,
                               fluid_6months_imp_demand_va_yrs),
                          ~extract_sl_continuous(..1, ..2, ..3, eye_demand, round = FALSE)) %>%
  set_names(c("va_yr_1", "va_yr_2", "va_yr_3")) %>%
  bind_rows(.id = "outcome")

