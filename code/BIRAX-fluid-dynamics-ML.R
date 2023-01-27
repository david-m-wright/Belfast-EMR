
source(rprojroot::find_rstudio_root_file("code/BIRAX-fluid-dynamics.R"))

# Setup dataset and learner stack

library(tlverse)
library(sl3)
library(conflicted)
library(pROC)
library(caret)
# library(fastshap)
library(reticulate)
library(polspline)
library(ranger)
library(xgboost)
library(caret)
library(Rsolnp)

# Increase java memory limit before calling bartMachine
options(java.parameters = "-Xmx30g") # 30Gb
library(bartMachine)
# use_condaenv("shap")


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


# Imputed data used for model fitting (both superlearner and conventional methods)
fluid_baseline_imp

# Multinomial regression for cluster
multi_fit <- nnet::multinom(as.formula(paste0("cluster ~ `", paste(noa_volumes_selected, collapse = "` + `"), "`")), data = fluid_baseline_imp) 

# Confusion matrix
multi_confusion <- confusionMatrix(as.factor(predict(multi_fit)), reference = as.factor(fluid_baseline_imp$cluster), dnn = c("Predicted Cluster", "True Cluster"), mode = "everything")

# Individual multinomial predictions for each variable
single_var_fit <- noa_volumes_selected %>% enframe(name = NULL) %>%
  mutate(fit_obj = map(
    .x = value,
    .f = ~ nnet::multinom(as.formula(paste0("cluster ~ `", .x, "`")), data = fluid_baseline_imp)
  )) %>%
  
  # Multiclass ROCs
  mutate(multi_roc = map(
    .x = fit_obj,
    .f = ~ multiclass.roc(fluid_baseline_imp$cluster, predict(.x, type = "prob"))
  )) %>% 
  
  # Multiclass AUCs (no CIs as these are a mean across multiple ROC curves)
  mutate(multi_auc = map_dbl(.x = multi_roc, .f = ~ .$auc))


## Superlearner models ##

# Create sl3 task
# Note imputation of missing values (median for continuous, mode for categorical) is
# performed in task setup if is has not been done already
cluster_task <- make_sl3_Task(
  data = fluid_baseline_imp,
  covariates = noa_volumes_selected,
  outcome = "cluster",
  outcome_type = "categorical",
  folds = 5
)


# cluster_task <- make_sl3_Task(
#   data = fluid_6months_imp,
#   covariates = paste(noa_volumes_selected, rep(c(0, 6), each = length(noa_volumes_selected)), sep = "_"),
#   outcome = "cluster",
#   outcome_type = "categorical",
#   folds = 5
# )
# 
# cluster_task <- make_sl3_Task(
#   data = fluid_6_12_imp,
#   covariates = paste(noa_volumes_selected, rep(c(0, 6, 12), each = length(noa_volumes_selected)), sep = "_"),
#   outcome = "cluster",
#   outcome_type = "categorical",
#   folds = 5
# )


# Learners that can work for this task
sl3_list_learners("categorical")
# sl_learner_descriptions <- read_csv(find_rstudio_root_file("Machine learning", "sl3_learner_descriptions.csv"))
# sl_variant_descriptions <- read_csv(find_rstudio_root_file("Machine learning", "sl3_variant_descriptions.csv"))

# Note that glm does not support multinomial outcome so cannot be used in its raw form
# Can use this method instead
# multinom_gf = make_learner(Lrnr_independent_binomial, make_learner(Lrnr_glm_fast)),

# Choose base learners (instantiate with R6 method $new())
lrn_mean <- Lrnr_mean$new()
lrn_lasso <- Lrnr_glmnet$new()
lrn_ridge <- Lrnr_glmnet$new(alpha = 0)
lrn_polspline <- Lrnr_polspline$new()
lrn_ranger100 <- Lrnr_ranger$new(num.trees = 100)
lrn_xgb_fast <- Lrnr_xgboost$new()
lrn_xgb50 <- Lrnr_xgboost$new(nrounds = 50)
lrn_nnet <- Lrnr_nnet$new()
lrn_bartMachine <- Lrnr_bartMachine$new(serialize =TRUE)
# Caret functions do not work
# lrn_caret_nnet <- Lrnr_caret$new(algorithm = "nnet")
# lrn_caret_bartMachine <- Lrnr_caret$new(algorithm = "bartMachine", method = "boot",  tuneLength = 10)


# For multinomial predictions are class probabilities

# Stack the learners so they are fitted simultaneously
learners <- c(lrn_mean, 
              lrn_lasso, 
              lrn_ridge, 
              lrn_polspline, 
              lrn_ranger100, 
              lrn_xgb_fast, 
              lrn_xgb50,
              lrn_nnet,
              lrn_bartMachine
              # lrn_caret_nnet,
              # lrn_caret_bartMachine
)

names(learners) <- c("mean", 
                     "lasso", 
                     "ridge", 
                     "polspline", 
                     "ranger_100", 
                     "xgb_fast", 
                     "xgb_50",
                     "nnet",
                     "bartMachine"
                     # "caret_nnet",
                     # "caret_bartMachine"
)



sl_learners <- learners[c(1:3, 5:7)]
#sl_learners <- learners

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
tictoc::tic()
sl_fit <- sl$train(cluster_task)
tictoc::toc()


################
# Function extract predictions from multinomial superlearner fit
# Args:
# sl_model = trained superlearner model
# sl_task = task on which the model was trained on
# id_cols = ID columns in the data input to the task
ExtractMultinomialPredictions <-
  function(sl_model, sl_task, id_cols) {
    # Predict probabilities of each class membership
    sl_preds_raw <- sl_model$predict() %>%
      map_dfr(unlist)
    
    outcome <- cluster_task$Y
    # Class with highest predicted probability
    pred_class <-
      # factor(
      apply(
        sl_preds_raw,
        1,
        FUN = function(x) {
          names(which.max(x))
        }
      )
    # , levels = levels(outcome))
    
    # Name output columns
    outcome_name <- cluster_task$nodes$outcome
    pred_outcome <- paste0("pred_", outcome_name)
    names(sl_preds_raw) <-
      paste("prob", outcome_name, names(sl_preds_raw), sep = "_")
    
    tibble(
      # ID columns
      cluster_task$get_data(columns = id_cols),
      # Outcome variable
      !!outcome_name := outcome,
      # Predicted class
      !!pred_outcome := pred_class,
      # Predicted class probabilities
      sl_preds_raw
    )
    
  }
# Still need to harmonise class of predicted classes

sl_out <- ExtractMultinomialPredictions(sl_fit, cluster_task, c("PatientID", "EyeCode")) %>% 
  mutate(across(c(cluster, pred_cluster), as.factor)) #%>% 

# Attempting to adjust cutoffs is not effective
# sl_out %>% 
# pivot_longer(cols = matches("prob")) %>% 
# group_by(PatientID, EyeCode, cluster, pred_cluster) %>%
# mutate(clus_rank = dense_rank(desc(value))) %>% 
# summarise(pred_cluster_second = as.factor(if_else(value[name == "prob_cluster_1"] > 0.6, "1", str_extract(name, "[0-9]")[clus_rank == 2])), .groups = "drop")
# confusionMatrix(sl_out$pred_cluster_second, reference = sl_out$cluster, dnn = c("Predicted cluster", "True cluster"), mode = "everything")

# Generate confusion matrix using predicted classes
sl_confusion <- confusionMatrix(sl_out$pred_cluster, reference = sl_out$cluster, dnn = c("Predicted cluster", "True cluster"), mode = "everything")
sl_confusion


sl_out %>% 
  pivot_longer(cols=matches("prob_")) %>% 
  filter(cluster == str_extract(name, "[0-9]")) %>% 
  ggplot(aes(y = value, x = cluster)) +
  geom_boxplot()


sl_out %>% 
  filter(cluster == 3) %>% 
  arrange(desc(prob_cluster_3))

