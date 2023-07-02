# Superlearner models to predict injection demand by treatment year

source(rprojroot::find_rstudio_root_file("code/BIRAX-injection-demand.R"))
set.seed(1983)

## Superlearner models ##

# 
# # Full multivariate model - not clear what this is estimating
# demand_task_mv_6 <- make_sl3_Task(
#   data = fluid_6months_imp_demand,
#   covariates = paste(noa_volumes_selected, rep(c(0, 6), each = length(noa_volumes_selected)), sep = "_"),
#   outcome = c("in_yr_1", "in_yr_2", "in_yr_3"),
#   outcome_type = "multivariate",
#   id = "PatientID",
#   folds = n_folds
# )
# 
# lrn_mean_mv <- Lrnr_multivariate$new(lrn_mean)
# lrn_glm_mv <- Lrnr_multivariate$new(lrn_glm)
# 
# stack_mv <- Stack$new(c(lrn_mean_mv, lrn_glm_mv))
# sl_mv <- Lrnr_sl$new(stack_mv)
# sl_fit_mv_6 <- sl_mv$train(demand_task_mv_6)



# Individual models for each year

# Create sl3 tasks for each year of injections
demand_task_yr1_6 <- make_sl3_Task(
  data = fluid_6months_imp_demand,
  covariates =  noa_covariates_0_6,
  outcome = "in_yr_1",
  outcome_type = "continuous",
  id = "PatientID",
  folds = n_folds
)

demand_task_yr2_6 <- make_sl3_Task(
  data = fluid_6months_imp_demand,
  covariates = noa_covariates_0_6,
  outcome = "in_yr_2",
  outcome_type = "continuous",
  id = "PatientID",
  folds = n_folds
)

demand_task_yr3_6 <- make_sl3_Task(
  data = fluid_6months_imp_demand,
  covariates = noa_covariates_0_6,
  outcome = "in_yr_3",
  outcome_type = "continuous",
  id = "PatientID",
  folds = n_folds
)

# List the demand tasks
demand_tasks_yrs <- list(demand_task_yr1_6, demand_task_yr2_6, demand_task_yr3_6) 
  

# Fit the superlearner
sl_fit_yrs_6 <- map(demand_tasks_yrs, ~ sl$train(.))

# Calculate the cross validated risk
sl_fit_yrs_6_cv_risk <-
  map(sl_fit_yrs_6, ~ .$cv_risk(eval_fun = loss_squared_error)) %>% 
  set_names(c("in_yr_1", "in_yr_2", "in_yr_3")) 

# Extract predictions
sl_preds_yrs_6  <-
  map2(
    demand_tasks_yrs,
    sl_fit_yrs_6,
    ~ extract_sl_continuous(.x, .y, fluid_6months_imp_demand, eye_demand, round = TRUE)
  ) %>%
  set_names(c("in_yr_1", "in_yr_2", "in_yr_3")) %>%
  bind_rows(.id = "outcome")


