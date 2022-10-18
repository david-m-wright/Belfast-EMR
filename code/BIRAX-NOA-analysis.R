# BIRAX project - analysis of NOA data

library(tidyverse)
library(rprojroot)
library(pROC)
library(caret)


source(find_rstudio_root_file("code/BIRAX-assemble-NOA-data.R"))


# Diagnostic test accuracy for each NOA variable #

# Binary variables - accuracy, sensitivity and specificity
# P. value is comparing the accuracy to the no information rate

# Get test statistics from confusion matrix
noa_conf_matrix <-  oct_visits %>%
  select(noa_fluid, noa_retinal) %>%
  select(where(is.factor)) %>%
  map( ~ confusionMatrix(.x,
                         reference = as.factor(as.numeric(
                           oct_visits$injected
                         )),
                         mode = "everything")) %>%
  map_dfr(broom::tidy, .id = "Variable") %>%
  filter(term %in% c("accuracy", "sensitivity", "specificity"))

# Output table
noa_diag_stats <- inner_join(
  noa_conf_matrix %>%
    select(Variable, term, estimate) %>%
    pivot_wider(names_from = term, values_from = estimate),
  
  noa_conf_matrix %>%
    filter(term == "accuracy") %>%
    select(Variable, p.value),
  by = c("Variable")
) %>%
  
  transmute(
    Variable,
    accuracy,
    P_accuracy = format.pval(p.value, eps = 0.001, digits = 2),
    sensitivity,
    specificity
  )




# To do
# 
# Find the optimum variable and threshold to ensure a given specificity is given.
# 
# Problem should be detecting scans where NO injection will be required (specificity >> 1?).
# 
# Partial ROC?
#   
#   Query for Notal? Are all scans centred on the same place?.
# 
# Categorical variables - confusion matrices, sensitivity and specificity
# 
# Multi-level model for thickness or volume.
# 
# Volume vs height
# 
# Secondary outcomes - does fluid volume correlate with VA?
#   Would there be a benefit to extending the analysis to before the first injection? i.e. to look at whether OCT indicates possible treatment outcomes later.
# Secondary analysis - how does NOA segmentation correlate with Heidelberg segmentation for each of the ETDRS regions?
#   To check - would there have been imaging only appointments or would all contain the possibility of injection on the same day?
#   