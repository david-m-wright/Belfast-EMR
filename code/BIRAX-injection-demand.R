# BIRAX - predicting number of injections required
library(rprojroot)
library(tidyverse)
library(janitor)
library(tlverse)
library(sl3)
library(conflicted)
library(pROC)
library(caret)
library(reticulate)
library(polspline)
library(ranger)
library(xgboost)
library(caret)
library(Rsolnp)
library(fastshap)
library(tictoc)


# Increase java memory limit before calling bartMachine
# options(java.parameters = "-Xmx30g") # 30Gb
# library(bartMachine)

library(reticulate)
use_condaenv("shap")


conflict_prefer("Stack", "sl3")
conflict_prefer("importance", "sl3")
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("intersect", "dplyr")
conflict_prefer("setdiff", "dplyr")
conflict_prefer("explain", "fastshap")
conflict_prefer("kable", "knitr")
conflict_prefer("slice", "dplyr")

# Setup parallel processing for training tasks
library(future)
ncores <- availableCores()-1
plan(multisession, workers = ncores)


# Load dataset
source(rprojroot::find_rstudio_root_file("code/BIRAX-fluid-dynamics.R"))

# Eyes to include, three years observation and > 3 fluid measurements
# This is a subset of the fluid_eyes list so can use imputed baseline fluid as is
eye_demand_raw <- eye %>% 
  left_join(three_yr_eyes %>% 
              mutate(exclude_three_yrs_observed = FALSE), by = c("PatientID", "EyeCode")) %>% 
  left_join(fluid_eyes %>% 
              select(PatientID, EyeCode) %>% 
              mutate(exclude_fluid_eyes = FALSE), by = c("PatientID", "EyeCode")) %>% 
  # Those with a baseline and a six month fluid measurement 
  left_join(fluid_6months_imp %>% 
               select(PatientID, EyeCode) %>% 
               mutate(exclude_6months_fluid = FALSE), by = c("PatientID", "EyeCode")) %>% 
  mutate(across(c(exclude_three_yrs_observed, exclude_fluid_eyes, exclude_6months_fluid), ~if_else(is.na(.), TRUE, .)))
eye_demand_raw %>% 
  count(across(matches("exclude")))

# eye_demand <- eye %>% 
#   inner_join(three_yr_eyes, by = c("PatientID", "EyeCode")) %>% 
#   inner_join(fluid_eyes, by = c("PatientID", "EyeCode"))

eye_demand <- eye_demand_raw %>% 
  filter(!exclude_three_yrs_observed, !exclude_fluid_eyes, !exclude_6months_fluid)

# Outcome - number of injections at 3 years
three_years_n_injections <- injections %>% 
  filter(years_treated <= 3) %>% 
  count(PatientID, EyeCode, name = "three_yr_injections")

# Baseline fluid measurements with outcome added
fluid_baseline_imp_demand <- fluid_baseline_imp %>% 
  inner_join(eye_demand, by = c("PatientID", "EyeCode")) %>% 
  inner_join(three_years_n_injections, by = c("PatientID", "EyeCode"))

fluid_baseline_imp_demand %>% 
  mutate(across(three_yr_injections, as.factor)) %>% 
  ggplot(aes(y = IRFVolumeNl, x = three_yr_injections)) +
  geom_violin() +
  lims(y = c(0, 0.5))


# Baseline and 6 month fluid measurements with outcome added
fluid_6months_imp_demand <- fluid_6months_imp %>% 
  inner_join(eye_demand, by = c("PatientID", "EyeCode")) %>% 
  inner_join(three_years_n_injections, by = c("PatientID", "EyeCode"))


# ROC performance of IRF volume alone
# pROC::multiclass.roc(response = fluid_baseline_imp_demand$three_yr_injections, predictor = fluid_baseline_imp_demand$IRFVolumeNl) 


# Subset fluid eye and VA tables for just those eyes in the demand analysis

fluid_history_demand <- fluid_history %>% 
  inner_join(eye_demand %>% 
               select(PatientID, EyeCode),
             by = c("PatientID", "EyeCode")) %>% 
  filter(months_since_index <= 36)

va_history_demand <- va_history %>% 
  inner_join(eye_demand %>% 
               select(PatientID, EyeCode),
             by = c("PatientID", "EyeCode")) %>% 
  filter(months_since_index <= 36)

thickness_history_demand <- thickness_history %>% 
  inner_join(eye_demand %>% 
               select(PatientID, EyeCode),
             by = c("PatientID", "EyeCode")) %>% 
  filter(months_since_index <= 36)


# injections_history_demand <- injections %>% 
#   filter(months_since_index <= 36)


## Superlearner models ##

# Create sl3 task
# Note imputation of missing values (median for continuous, mode for categorical) is
# performed in task setup if is has not been done already
demand_task <- make_sl3_Task(
  data = fluid_baseline_imp_demand,
  covariates = noa_volumes_selected,
  outcome = "three_yr_injections",
  outcome_type = "continuous",
  id = "PatientID",
  folds = 5
)


demand_task_6 <- make_sl3_Task(
  data = fluid_6months_imp_demand,
  covariates = paste(noa_volumes_selected, rep(c(0, 6), each = length(noa_volumes_selected)), sep = "_"),
  outcome = "three_yr_injections",
  outcome_type = "continuous",
  id = "PatientID",
  folds = 5
)
# 
# cluster_task <- make_sl3_Task(
#   data = fluid_6_12_imp,
#   covariates = paste(noa_volumes_selected, rep(c(0, 6, 12), each = length(noa_volumes_selected)), sep = "_"),
#   outcome = "cluster",
#   outcome_type = "categorical",
#   folds = 5
# )




# Learners that can work for this task
sl3_list_learners("continuous") %>% enframe() %>% data.frame()
# sl_learner_descriptions <- read_csv(find_rstudio_root_file("Machine learning", "sl3_learner_descriptions.csv"))
# sl_variant_descriptions <- read_csv(find_rstudio_root_file("Machine learning", "sl3_variant_descriptions.csv"))

# Note that glm does not support multinomial outcome so cannot be used in its raw form
# Can use this method instead
# multinom_gf = make_learner(Lrnr_independent_binomial, make_learner(Lrnr_glm_fast)),

# Choose base learners (instantiate with R6 method $new())
lrn_mean <- Lrnr_mean$new()
lrn_glm <- Lrnr_glm$new()
lrn_lasso <- Lrnr_glmnet$new()
lrn_ridge <- Lrnr_glmnet$new(alpha = 0)
lrn_polspline <- Lrnr_polspline$new()
lrn_ranger1000 <- Lrnr_ranger$new(num.trees = 1000)
lrn_xgb_fast <- Lrnr_xgboost$new()
lrn_xgb100 <- Lrnr_xgboost$new(nrounds = 100)
lrn_nnet <- Lrnr_nnet$new()
lrn_bartMachine <- Lrnr_bartMachine$new(serialize =TRUE)
lrn_hal9001 <- Lrnr_hal9001$new()
# Caret functions do not work
# lrn_caret_nnet <- Lrnr_caret$new(algorithm = "nnet")
# lrn_caret_bartMachine <- Lrnr_caret$new(algorithm = "bartMachine", method = "boot",  tuneLength = 10)
lrn_nnet_autotune <- Lrnr_caret$new(algorithm = "nnet", name = "NNET_autotune")

# For multinomial predictions are class probabilities

# Stack the learners so they are fitted simultaneously
learners <- c(lrn_mean, 
              lrn_glm,
              lrn_lasso, 
              lrn_ridge, 
              lrn_polspline, 
              lrn_ranger1000, 
              lrn_xgb_fast, 
              lrn_xgb100,
              lrn_nnet,
              # lrn_bartMachine,
               # lrn_hal9001
              # lrn_caret_nnet,
              # lrn_caret_bartMachine
              lrn_nnet_autotune
)

names(learners) <- c("mean", 
                     "glm",
                     "lasso", 
                     "ridge", 
                     "polspline", 
                     "ranger_1000", 
                     "xgb_fast", 
                     "xgb_100",
                     "nnet",
                     # "bartMachine",
                     # "hal9001"
                     # "caret_nnet",
                     # "caret_bartMachine"
                     "caret_nnet_autotune"
)



# sl_learners <- learners[c(1:3, 5:7)]
sl_learners <- learners

sl_learners %>% 
  enframe() %>% 
  transmute(Algorithm = name)
 
# Make the Super Learner 

# Learner stack should be trained on the full data (no cross validation)
# Stack the learners
# stack <- make_learner_stack(sl_learners)
stack <- Stack$new(sl_learners)


# # New syntax for stacking learners
# stack <- Stack$new(
#   "Lrnr_mean", 
#   "Lrnr_glmnet",
#   list("Lrnr_glmnet", 
#        alpha = 0),
#   list("Lrnr_ranger",
#        num.trees = 100),
#   "Lrnr_xgboost",
#   list("Lrnr_xgboost", 
#        nrounds = 50)
# )


# Meta learner should be trained on the cross validated predictions
# Use default metalearner (solnp)
sl <- Lrnr_sl$new(stack)

#  Train the superlearner
# 2 hrs
# tictoc::tic()
# sl_fit <- sl$train(demand_task)
# tictoc::toc()

# 4 hrs
tictoc::tic()
sl_fit_6 <- sl$train(demand_task_6)
tictoc::toc()


# sl_fit_cv_risk <- sl_fit$cv_risk()
sl_fit_6_cv_risk <- sl_fit_6$cv_risk(eval_fun = loss_squared_error)

# save(sl_fit, file = find_rstudio_root_file("workspaces/sl-fit.RData"))
# save(sl_fit_6, file = find_rstudio_root_file("workspaces/sl-fit-6.RData"))

# # Extract predictions 
# # Note that fractional numbers of injections can be predicted so round to nearest integer
# sl_preds <- tibble(demand_task$data,
#   pred_three_yr_injections = round(sl_fit$predict()))
# 
# # Strong positive association between predicted and true values
# cor.test(sl_preds$three_yr_injections, sl_preds$pred_three_yr_injections)
# summary(lm(three_yr_injections~pred_three_yr_injections, data = sl_preds))
# xtabs(data = sl_preds, ~pred_three_yr_injections+ three_yr_injections)


# Extract predictions 
# Note that fractional numbers of injections can be predicted so round to nearest integer
sl_preds_6 <- tibble(demand_task_6$data,
                   pred_three_yr_injections = round(sl_fit_6$predict()),
                   deviation = three_yr_injections - pred_three_yr_injections,
                   abs_deviation = abs(deviation)) %>% 
  bind_cols(eye_demand %>% 
              select(-PatientID)) %>% 
  mutate(
                   # Positional matching to original data
                year_start = as.factor(data.table::year(index_date)))
  

tabyl(sl_preds_6, three_yr_injections, pred_three_yr_injections)

# Strong positive association between predicted and true values
cor_sl_preds_6 <- cor.test(sl_preds_6$three_yr_injections, sl_preds_6$pred_three_yr_injections)

# Calibration slope
cal_6 <- lm(three_yr_injections~pred_three_yr_injections, data = sl_preds_6)
summary(cal_6)


# save(sl_preds, file = find_rstudio_root_file("workspaces/sl-preds.RData"))
# save(sl_preds_6, file = find_rstudio_root_file("workspaces/sl-preds-6.RData"))


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
    folds = 5
  )
  # Predict number of injections
  learner_fit_predict(object, task = sl_wrap) 
}
# Test the prediction wrapper
# demand_prediction_wrapper(sl_fit_6, fluid_6months_imp_demand[, c("PatientID", paste(noa_volumes_selected, rep(c(0, 6), each = length(noa_volumes_selected)), sep = "_")), with =FALSE])

# Calculate SHAP values
# 30 mins for 10 simulations
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
    nsim = 10,
    adjust = TRUE
  ) 
toc()

