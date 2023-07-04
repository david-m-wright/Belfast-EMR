# Validation of injection demand model using TASMC data

# Setup environment for superlearner
library(rprojroot)
library(tidyverse)
library(data.table)
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
options(java.parameters = "-Xmx30g") # 30Gb
# options(java.parameters = "-Xmx1g") # 1Gb
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

# Assemble the analysis cohort
source(find_rstudio_root_file("Code/TASMC-assemble-cohort.R"))

# Load the trained prediction models and the scaling coefficients

# Imputation of missing variables
noa_fluid
noa_general

# All volume measurements except total volumes, as these are linear combinations of IRF and SRF
noa_volumes <- noa %>% 
  select(matches("volume|_vol"), -matches("TRF")) %>% 
  names()

# Quantify missingness for volumes
fluid_history %>% 
  select(all_of(noa_volumes)) %>% 
  summarise(across(.cols = all_of(noa_volumes), .fns = ~sum(is.na(.))/length(.)*100)) %>% 
  pivot_longer(cols = everything()) %>% 
  arrange(desc(value)) %>% print(n=Inf)

# Just select volumes with a negligible amount of imputation
# noa_volumes_selected <- noa %>% 
#   select(all_of(noa_volumes), -matches("Superior6|Inferior6")) %>% 
#   names()
# Or include them all
noa_volumes_selected <- noa_volumes

# Impute missing volumes with zero 
# Scale by mean and sd from baseline measurements in the Belfast dataset
imputation_scaling_values <- read_rds(find_rstudio_root_file("SL-models", "imputation-scaling-values.RDS"))

fluid_history_imp <- eye %>% 
  select(PatientID, EyeCode, eye_id) %>% 
  inner_join(
  bind_cols(
  fluid_history %>%
    select(-all_of(noa_volumes)),
  
  fluid_history %>%
    select(all_of(noa_volumes)) %>%
    apply(
      2,
      FUN = function(x)
        if_else(is.na(x), 0, x)
    ) %>%
    sweep(2, STATS = imputation_scaling_values$mean_training, FUN = "-") %>%
    sweep(2, STATS = imputation_scaling_values$sd_training, FUN = "/")
),
by = c("PatientID", "EyeCode")) %>%
  
  arrange(PatientID, EyeCode, months_since_index) %>%
  group_by(PatientID, EyeCode) %>%
  mutate(oct_series = row_number()) %>%
  ungroup() %>%
  as.data.table()

# Just those with a snapshot at 6 months
fluid_6months_imp <- fluid_history_imp %>% 
  filter(snapshot, follow_up_month <= 6) %>% 
  select(PatientID, EyeCode, eye_id, follow_up_month, all_of(noa_volumes_selected)) %>% 
  pivot_wider(id_cols = c(PatientID, EyeCode), 
              names_from = follow_up_month,
              values_from = all_of(noa_volumes_selected)) %>% 
  na.omit()

# Eyes to include:
# outcome recorded: three years observation? Not for 1 and 2 year models, just predict for fewer eyes
# and > 3 fluid measurements (short sequences unlikely to be informative) *** NOT USED
eye_demand_raw <- eye %>% 
  mutate(exclude_three_yrs_observed = years_observed >= 3) %>% 
  
  # Number of fluid measurements
  left_join(fluid_history %>% 
  count(PatientID, EyeCode, name = "measurements") %>% 
  filter(measurements > 3) %>% 
    mutate(exclude_fluid_eyes = FALSE)
  , by = c("PatientID", "EyeCode")) %>%
  
  # Those with a baseline and a six month fluid measurement 
  left_join(fluid_6months_imp %>% 
              select(PatientID, EyeCode) %>% 
              mutate(exclude_6months_fluid = FALSE), by = c("PatientID", "EyeCode")) %>% 
  mutate(across(c(exclude_three_yrs_observed, exclude_fluid_eyes, exclude_6months_fluid), ~if_else(is.na(.), TRUE, .)))

eye_demand_raw %>% 
  count(across(matches("exclude")))

eye_demand <- eye_demand_raw %>% 
  filter(
    # Relax the three year requirement? - excludes half of the TASMC patients
     # !exclude_three_yrs_observed, 
    # !exclude_fluid_eyes, 
    !exclude_6months_fluid)

# Outcomes - 
# Numbers of injections by treatment year
injections_yr <- injections %>% 
  count(PatientID, EyeCode, treatment_year) %>% 
  pivot_wider(names_from = treatment_year, values_from = n, names_prefix = "in_yr_", values_fill = 0) %>% 
  # number of injections at 3 years
  mutate(three_yr_injections = in_yr_1 + in_yr_2 + in_yr_3)


# Baseline and 6 month fluid measurements with outcome added
fluid_6months_imp_demand <- fluid_6months_imp %>% 
  inner_join(eye_demand, by = c("PatientID", "EyeCode")) %>% 
  # inner_join(three_years_n_injections, by = c("PatientID", "EyeCode"))
  inner_join(injections_yr, by = c("PatientID", "EyeCode")) 


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


# Fluid measurements to include as predictors in superlearner models
# Baseline and six month
noa_covariates_0_6 <- paste(noa_volumes_selected, rep(c(0, 6), each = length(noa_volumes_selected)), sep = "_")


# Helper functions for handling superlearner fits and predictions

# Extract predictions from superlearner fit and add eye code as another identifier
# Args:
# sl_task = sl3 task 
# sl_fit = trained sl3 model
# sl_dat = data input to sl3 task but also containing eye code
# eye_chars = data.table of additional eye characteristics (keys: PatientID, EyeCode)
# round = logical indicating whether to round the predictions before calculating the deviations (e.g. for integer predictions)
extract_sl_continuous <- function(sl_task, sl_fit, sl_dat, eye_chars, round = FALSE) {
  
  preds_Y <- if(round){
    round(sl_fit$predict(sl_task))
  } else {
    sl_fit$predict(sl_task)
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
    # bind_cols(sl_dat %>%
    #             select(EyeCode)) %>%
    # Additional characteristics
    inner_join(eye_chars, by = "eye_id") %>%
    mutate(year_start = as.factor(data.table::year(index_date)))
}


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