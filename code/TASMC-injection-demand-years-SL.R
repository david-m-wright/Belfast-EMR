# Superlearner models to predict injection demand at end of each treatment year

source(rprojroot::find_rstudio_root_file("code/TASMC-injection-demand.R"))
set.seed(1983)

# Individual models for each year

# Create sl3 tasks for each year of injections
demand_task_yr1_6 <- make_sl3_Task(
  # data = filter(fluid_6months_imp_demand, years_observed >=1),
  data = fluid_6months_imp_demand,
  covariates =  noa_covariates_0_6,
  outcome = "in_yr_1",
  outcome_type = "continuous",
  id = "eye_id"
)

demand_task_yr2_6 <- make_sl3_Task(
  # data = filter(fluid_6months_imp_demand, years_observed >= 2),
  data = fluid_6months_imp_demand,
  covariates = noa_covariates_0_6,
  outcome = "in_yr_2",
  outcome_type = "continuous",
  id = "eye_id"
)

demand_task_yr3_6 <- make_sl3_Task(
  # data = filter(fluid_6months_imp_demand, years_observed >= 3),
  data = fluid_6months_imp_demand,
  covariates = noa_covariates_0_6,
  outcome = "in_yr_3",
  outcome_type = "continuous",
  id = "eye_id"
)

# List the demand tasks
demand_tasks_yrs <- list(demand_task_yr1_6, demand_task_yr2_6, demand_task_yr3_6) 


# Load the fitted model
sl_fit_yrs_6 <- read_rds(find_rstudio_root_file("SL-models", "sl_fit_yrs_6.RDS"))

# Extract predictions
sl_preds_yrs_6  <-
  map2(
    demand_tasks_yrs,
    sl_fit_yrs_6,
    ~ extract_sl_continuous(.x, .y, fluid_6months_imp_demand, eye_demand, round = TRUE)
  ) %>%
  set_names(c("in_yr_1", "in_yr_2", "in_yr_3")) %>%
  bind_rows(.id = "outcome")

