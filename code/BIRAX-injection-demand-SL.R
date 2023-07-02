# Superlearner models to predict injection demand at 3 years

source(rprojroot::find_rstudio_root_file("code/BIRAX-injection-demand.R"))
set.seed(1983)

# Create sl3 task
# Note imputation of missing values (median for continuous, mode for categorical) is
# performed in task setup if is has not been done already
demand_task_6 <- make_sl3_Task(
  data = fluid_6months_imp_demand,
  # covariates = c("index_age", "Gender", "va_logmar", noa_covariates_0_6),
  covariates = noa_covariates_0_6,
  outcome = "three_yr_injections",
  outcome_type = "continuous",
  id = "PatientID",
  folds = n_folds
)


#  Train the superlearner and calculate cross validated risks
tic()
sl_fit_6 <- sl$train(demand_task_6)
sl_fit_6_cv_risk <- sl_fit_6$cv_risk(eval_fun = loss_squared_error)
toc()

# save(sl_fit, file = find_rstudio_root_file("workspaces/sl-fit.RData"))
# save(sl_fit_6, file = find_rstudio_root_file("workspaces/sl-fit-6.RData"))


sl_preds_6  <- extract_sl_continuous(demand_task_6, sl_fit_6, fluid_6months_imp_demand, eye_demand, round = TRUE)
# save(sl_preds_6, file = find_rstudio_root_file("workspaces/sl-preds-6.RData"))


# Function to summarise predictions from a superlearner fit 
# Args:
# sl_predictions = predictions extracted using extract_sl_continuous()
summarise_sl_fit <- function(sl_predictions){
  bind_rows(
  cor.test(sl_predictions$Y, sl_predictions$pred_Y) %>% 
    tidy(),  
  
  lm(Y~pred_Y, data = sl_predictions) %>% 
    tidy() %>% 
    mutate(method = "Linear regression (calibration)"),
  
  sl_predictions %>% summarise(Deviation = mean(deviation),
                               `Absolute deviation` = mean(abs_deviation)) %>% 
    pivot_longer(cols = everything(), names_to = "method", values_to = "estimate") %>% 
    mutate(term = "mean")
  
  ) %>% 
     transmute(method, term, estimate, p = format.pval(p.value, eps = 0.001))
    
}

summarise_sl_fit(sl_preds_6)


## Interpretability

## Setup prediction wrapper that will be called by explain()
# This produces a prediction for each of the Frankenstein matrices constructed by fastSHAP
# Can have only two arguments, object and newdata
demand_6_prediction_wrapper <- function(object, newdata){
  sl_wrap <- make_sl3_Task(
    data = cbind(three_yr_injections = NA, newdata),
    covariates = paste(noa_volumes_selected, rep(c(0, 6), each = length(noa_volumes_selected)), sep = "_"),
    outcome = "three_yr_injections",
    outcome_type = "continuous",
    id = "PatientID",
    folds = n_folds
  )
  # Predict number of injections
  learner_fit_predict(object, task = sl_wrap) 
}
# Test the prediction wrapper
# demand_prediction_wrapper(sl_fit_6, fluid_6months_imp_demand[, c("PatientID", paste(noa_volumes_selected, rep(c(0, 6), each = length(noa_volumes_selected)), sep = "_")), with =FALSE])

# Calculate SHAP values
# Setup parallel backend for SHAP values
# library(doParallel)
# cl <- makeCluster(ncores)
# registerDoParallel(cl)


# 30 mins for 10 simulations (series)
# 6 mins for 2
tic()
sl_shap_all <-
  explain(
    sl_fit_6,
    # The input matrix must be the dataset that is input to the sl3 task
    # rather than the one output by the task where imputation has been automatically carried out by sl3
    # Using the output from the sl3 task leads to predictions not matching fitted values and failure of local additivity of SHAP values
    #  X = as.data.frame(sl_dat[, sl_selected]),
    # X = as.data.frame(fluid_6months_imp_demand,[, sl_selected]),
    X = as.data.frame(fluid_6months_imp_demand[, c("PatientID", paste(noa_volumes_selected, rep(c(0, 6), each = length(noa_volumes_selected)), sep = "_")), with =FALSE]),
    # X = imp,
    pred_wrapper = demand_6_prediction_wrapper,
    nsim = 2,
    adjust = TRUE#,
    # .parallel = TRUE
  ) 
toc()