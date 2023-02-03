# Superlearner models to predict treatment duration

source(rprojroot::find_rstudio_root_file("code/BIRAX-injection-demand.R"))

# Define task
treatment_task_6 <- make_sl3_Task(
  data = fluid_6months_imp_demand,
  covariates =  noa_covariates_0_6,
  outcome = "years_treated",
  outcome_type = "continuous",
  id = "PatientID",
  folds = n_folds
)

# Train superlearner
tictoc::tic()
sl_fit_treatment_6 <- sl$train(treatment_task_6)
sl_fit_treatment_6_cv_risk <- sl_fit_treatment_6$cv_risk(eval_fun = loss_squared_error)
tictoc::toc()

# Extract predictions
sl_preds_treatment_6  <- extract_sl_continuous(treatment_task_6, sl_fit_treatment_6, fluid_6months_imp_demand, eye_demand, round = FALSE)


