# BIRAX - predicting number of injections required

# Setup environment for superlearner
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
library(broom)

# Increase java memory limit before calling bartMachine
# options(java.parameters = "-Xmx30g") # 30Gb
options(java.parameters = "-Xmx1g") # 1Gb
library(bartMachine)

# Call python environment with SHAP package for force plots of SHAP values
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

eye_demand <- eye_demand_raw %>% 
  filter(
    # Relax the three year requirement?
     !exclude_three_yrs_observed, 
         !exclude_fluid_eyes, 
         !exclude_6months_fluid)

# Outcomes - 
# Numbers of injections by treatment year
injections_yr <- injections %>% 
  count(PatientID, EyeCode, treatment_year) %>% 
  pivot_wider(names_from = treatment_year, values_from = n, names_prefix = "in_yr_", values_fill = 0) %>% 
# number of injections at 3 years
  mutate(three_yr_injections = in_yr_1 + in_yr_2 + in_yr_3)

# three_years_n_injections <- injections %>%
#   filter(years_treated <= 3) %>%
#   count(PatientID, EyeCode, name = "three_yr_injections")

  
# Baseline fluid measurements with outcome added
fluid_baseline_imp_demand <- fluid_baseline_imp %>% 
  inner_join(eye_demand, by = c("PatientID", "EyeCode")) %>% 
  # inner_join(three_years_n_injections, by = c("PatientID", "EyeCode"))
  inner_join(injections_yr, by = c("PatientID", "EyeCode"))

fluid_baseline_imp_demand %>% 
  mutate(across(three_yr_injections, as.factor)) %>% 
  ggplot(aes(y = IRFVolumeNl, x = three_yr_injections)) +
  geom_violin() +
  lims(y = c(0, 0.5))


# Baseline and 6 month fluid measurements with outcome added
fluid_6months_imp_demand <- fluid_6months_imp %>% 
  inner_join(eye_demand, by = c("PatientID", "EyeCode")) %>% 
  # inner_join(three_years_n_injections, by = c("PatientID", "EyeCode"))
  inner_join(injections_yr, by = c("PatientID", "EyeCode"))

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



# Find last VA recorded for each treatment year
va_yr <- va_history_demand[va_history_demand[,.I[years_since_index == max(years_since_index)], by=.(PatientID, EyeCode, treatment_year)]$V1][
  treatment_year != 0, .(PatientID, EyeCode, treatment_year, va_logmar)] %>% 
  pivot_wider(names_from = treatment_year, values_from = va_logmar, names_prefix = "va_logmar_yr_")

fluid_6months_imp_demand_va_yr1 <- fluid_6months_imp_demand %>% 
  inner_join(va_yr %>% 
               select(PatientID, EyeCode, va_logmar_yr_1) %>% 
               filter(!is.na(va_logmar_yr_1)), by = c("PatientID", "EyeCode"))

fluid_6months_imp_demand_va_yr2 <- fluid_6months_imp_demand %>% 
  inner_join(va_yr %>% 
               select(PatientID, EyeCode, va_logmar_yr_2) %>% 
               filter(!is.na(va_logmar_yr_2)), by = c("PatientID", "EyeCode"))

fluid_6months_imp_demand_va_yr3 <- fluid_6months_imp_demand %>% 
  inner_join(va_yr %>% 
               select(PatientID, EyeCode, va_logmar_yr_3) %>% 
               filter(!is.na(va_logmar_yr_3)), by = c("PatientID", "EyeCode"))

fluid_6months_imp_demand_va_yrs <- list(fluid_6months_imp_demand_va_yr1, fluid_6months_imp_demand_va_yr2, fluid_6months_imp_demand_va_yr3)

# va_history %>% 
#   filter(PatientID == "0159B375-C58B-B305-D2D3-BBF28C638000", EyeCode == "L")


# Fluid measurements to include as predictors in superlearner models
# Baseline and six month
noa_covariates_0_6 <- paste(noa_volumes_selected, rep(c(0, 6), each = length(noa_volumes_selected)), sep = "_")


## Define SuperLearner algorithm for predicting continuous variables

# Learners that can work for this task
sl3_list_learners("continuous") %>% enframe() %>% data.frame()
sl_learner_descriptions <- read_csv(find_rstudio_root_file("data-dictionary", "sl3_learner_descriptions.csv"))
sl_variant_descriptions <- read_csv(find_rstudio_root_file("data-dictionary", "sl3_variant_descriptions.csv"))

# Base learners (instantiate with R6 method $new())
lrn_mean <- Lrnr_mean$new()
lrn_glm <- Lrnr_glm$new()
lrn_lasso <- Lrnr_glmnet$new()
lrn_ridge <- Lrnr_glmnet$new(alpha = 0)
lrn_polspline <- Lrnr_polspline$new()
lrn_ranger1000 <- Lrnr_ranger$new(num.trees = 1000)
lrn_xgb_fast <- Lrnr_xgboost$new()
lrn_xgb100 <- Lrnr_xgboost$new(nrounds = 100)
lrn_nnet <- Lrnr_nnet$new()
lrn_nnet_autotune <- Lrnr_caret$new(algorithm = "nnet", name = "NNET_autotune")
lrn_bartMachine <- Lrnr_bartMachine$new(serialize =TRUE) # Must serialise to save trained models
lrn_hal9001 <- Lrnr_hal9001$new()

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
              lrn_nnet_autotune,
              lrn_bartMachine,
              lrn_hal9001
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
                     "caret_nnet_autotune",
                     "bartMachine",
                     "hal9001"
)


sl_learners <- learners[c(1:8, 10)]
sl_learners %>% 
  enframe() %>% 
  transmute(Algorithm = name)

# Stack the learners
stack <- Stack$new(sl_learners)

# Make the Super Learner using default meta-learner (solnp)
# Ensures meta learner is trained on the cross validated predictions
sl <- Lrnr_sl$new(stack)

# Number of folds to use for superlearner cross-validation
n_folds <- 5

# Extract predictions from superlearner fit and add eye code as another identifier
# Args:
# sl_task = sl3 task 
# sl_fit = trained sl3 model
# sl_dat = data input to sl3 task but also containing eye code
# eye_chars = data.table of additional eye characteristics (keys: PatientID, EyeCode)
# round = logical indicating whether to round the predictions before calculating the deviations (e.g. for integer predictions)
extract_sl_continuous <- function(sl_task, sl_fit, sl_dat, eye_chars, round = FALSE) {
  
  preds_Y <- if(round){
    round(sl_fit$predict())
  } else {
    sl_fit$predict()
  }
    
  tibble(
    sl_task$data,
    Y = sl_task$Y,
    pred_Y = preds_Y,
    # Calculate deviations between prediction and observed value
    deviation = Y - pred_Y,
    abs_deviation = abs(deviation)
  ) %>%
    # Positional matching to original data to add secondary ID (EyeCode)
    bind_cols(sl_dat %>%
                select(EyeCode)) %>%
    # Additional characteristics
    inner_join(eye_chars, by = c("PatientID", "EyeCode")) %>%
    mutate(year_start = as.factor(data.table::year(index_date)))
}


