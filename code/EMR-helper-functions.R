# Helper functions
library(janitor)
library(tidyverse)

# Function to clean up NISA study IDs
# Args:
# x = character vector of IDs to be cleaned
# Value: character vector of cleaned IDs
CleanNISA <- function(x){
  toupper(str_replace_all(x, " ", "")) %>% 
    # For NICOLA IDs
    if_else(str_detect(., "^N[0-9]{5}"), .,
            # For NISA IDs
            str_replace(., "(2)([0-9]{3})", "\\1-\\2") %>% 
              str_extract("^(NISA|DMO|NISA2-|GWN)[0-9]{3}"))
}


# Function to select most frequently occuring non-missing value from a set of values
# Or first (using sort()) in case of ties
# Args:
# x = vector of values to count votes in
# Value: value with the most votes
MaxVotes <- function(x){
  if(all(is.na(x))){return(NA)}
  names(sort(table(x), decreasing = T))[1]
}


# Functions for generating descriptive statistics
# Continous variables


# Function to summarise continuous variables
SummariseContinuous <- function(x){
  tibble(Min = min(x, na.rm = T),
         Median = median(x, na.rm = T),
         Mean = mean(x, na.rm = T),
         Max = max(x, na.rm = T),
         SD = sd(x, na.rm = T),
         Missing = sum(is.na(x))) %>% 
    mutate_all(round, digits = 2)
}


# SD and mean for continous variables
SummariseSD <- function(x, digits = 1, n = FALSE){
  out <- paste0(round(mean(x, na.rm = T), digits), 
         " (", 
         round(sd(x, na.rm = T), digits), 
         ")")
  # Report number of valid records
  if(n){
   out <- paste0(out, " [", 
          sum(!is.na(x)),
          "]")
  }
  out
}

# Median and interquartile range for continuous variables
SummariseMedian <- function(x, digits = 1, n = FALSE){
  out <- paste0(round(median(x, na.rm = T), digits), 
         " (", 
         round(quantile(x, probs = c(0.25), na.rm = T), digits), 
         ", ",
         round(quantile(x, probs = c(0.75), na.rm = T), digits), 
         ")")
  if(n){
    out <- paste0(out, " [", 
           sum(!is.na(x)),
           "]")
  }
  out    
}

# Counts and percentages for discrete variables
# Args: dat = data frame from which to tabulate variables
# var1_list = vars() specification listing variables to tabulate (rows)
SummariseDiscreteOneway <- function(dat, var1_list){
  
  map_dfr(var1_list, ~tabyl(dat = {{dat}}, var1 = !!.) %>% 
            adorn_pct_formatting(digits = 1) %>%
            
            # Add columns for variable and level names
            mutate(Variable = names(.)[1]) %>% 
            rename(Level = names(.)[1])
  ) %>% 
    mutate(Value = paste0(n, " (", percent, ")")) %>% 
    select(Variable, Level, Value)
}

# patients %>% 
#   SummariseDiscreteOneway(vars(Gender, EthnicCode))





# Twoway display of continuous variables
SummariseContinousTwoway <- function(dat, var1_list, var2, summary_function = SummariseSD){
  # Which type of summary to take 
  summary_function <- summary_function
  
  # Select the data to summarise
  selected_data <- if("quosures" %in% class(var1_list)){
    {{dat}} %>% 
      select(!!!var1_list, {{var2}})
  } else {
    {{dat}} %>% 
      select(!!var1_list, {{var2}})
  }
  selected_data %>% 
    group_by({{var2}}) %>% 
    summarise_at({{var1_list}}, summary_function) %>%
    pivot_longer(-{{var2}}, names_to = "Variable") %>% 
    pivot_wider(names_from = {{var2}}, values_from = value)
}
# SummariseContinousTwoway(study_eye_data, c("Age", "HbA1c"), DM_DR_status)
# SummariseContinousTwoway(study_eye_data, c("Age", "HbA1c"), DM_DR_status, SummariseMedian)

# Perform t-test for each of a set of continuous variables,
# classified by a single discrete variable
# Args: dat = data frame from which to test variables
# var1_list = character vector or vars() specificiation listing variables to test (rows)
# var2 = unquoted name of variable to test by (columns)
SummariseTTest <- function(dat, var1_list, var2){
  
  vars_to_test <- if("quosures" %in% class(var1_list)){
    var1_list %>% 
      map(., quo_name)
  } else {
    var1_list
  }
  
  var2 <- enquo(var2)
  
  # Perform the t-tests
  vars_to_test %>% 
    map_dfr(~broom::tidy(t.test(as.formula(paste(., "~", quo_name(var2))), data = {{dat}}))) %>% 
    transmute(Variable = vars_to_test, 
              `t statistic` = round(statistic, digits = 2),
              df = round(parameter, digits = 1),
              P = format.pval(p.value, eps = 0.001, digits = 2))
}



# Discrete variables

# Frequencies and percentages for discrete variables
SummariseDiscrete <- function(x){
  as.vector(x) %>% 
    tabyl() %>% 
    rename(Value = ".") %>% 
    adorn_pct_formatting()
}


# Twoway frequency and percentage table for discrete variables
# Args: dat = data frame from which to tabulate variables
# var1_list = character vector or vars() specification listing variables to tabulate (rows)
# var2 = unquoted name of variable to tabulate by (columns)
SummariseDiscreteTwoway <- function(dat, var1_list, var2){
  
  # Generate the base frequency tables, one for each variable in the var1_list
  # Tabulates frequencies of var1 by var2
  cross_tabs <- if("quosures" %in% class(var1_list)){
    {{var1_list}}
  } else {
    var1_list %>% 
      syms() 
  }
  
  map_dfr(cross_tabs, ~tabyl(dat = {{dat}}, var1 = !!.,  var2 = {{var2}}) %>% 
            
            # Format percentages
            adorn_totals() %>%
            adorn_percentages(denominator = "col") %>%
            adorn_pct_formatting(digits = 1) %>%
            adorn_ns(position = "front") %>%
            
            # Add columns for variable and level names
            mutate(Variable = names(.)[1]) %>% 
            rename(Level = names(.)[1])) %>%
    
    # Drop the excess totals rows
    mutate(row_number = row_number()) %>% 
    filter(is.na(Level) | !(Level == "Total" & row_number != max(row_number))) %>% 
    mutate(Variable = if_else(row_number == max(row_number), "Total", Variable),
           Level = if_else(row_number == max(row_number), "-", Level)) %>% 
    select(-row_number) %>% 
    select(Variable, Level, everything())
  
}
#SummariseDiscreteTwoway(choroid_eyes, vars(Sex, Study_eye), DM_DR_status)
#SummariseDiscreteTwoway(choroid_eyes, c("Sex", "Study_eye"), DM_DR_status)



# Perform chi-squared test for each of a set of discrete variables,
# classified by a single discrete variable
# Args: dat = data frame from which to test variables
# var1_list = character vector or vars() specification listing variables to test (rows)
# var2 = unquoted name of variable to test by (columns)
SummariseChisqTest <- function(dat, var1_list, var2){
  
  cross_tabs <- if("quosures" %in% class(var1_list)){
    {{var1_list}}
  } else {
    var1_list %>% 
      syms() 
  }
  
  # Generate row names
  bind_cols(Variable = map_chr(cross_tabs, quo_name), 
            
            # Generate the base frequency tables, one for each variable in the var1_list
            # Tabulates frequencies of var1 by var2
            map_dfr(cross_tabs, ~tabyl(dat = {{dat}}, var1 = !!.,  var2 = {{var2}},
                                       show_missing_levels = FALSE,
                                       show_na = FALSE) %>% 
                      
                      # Peform the chi-squared tests
                      janitor::chisq.test() %>%
                      broom::tidy() %>%
                      
                      transmute(`Chi-sq statistic` = round(statistic, digits = 2),
                                df = round(parameter, digits = 1),
                                P = format.pval(p.value, eps = 0.001, digits = 2))) 
  )
  
}
# SummariseChisqTest(mtcars, c("cyl", "carb"), gear)
SummariseChisqTest(mtcars, vars(cyl, carb), gear)




# Function to display the diagnostic peformance stats of a predictor variable
# Args:
# dat = tibble containing data
# response = column containing response variable
# predictor = column containing predictor variable
# Value: tibble containing the summary stats
# Note that when this is applied to a grouped tibble the grouping variable is added automatically
# (with a message)
DiagnosticStats <- function(dat, response, predictor){
  response = enquo(response)
  predictor = enquo(predictor)
  
  dat %>%
    count(!!response, !!predictor) %>%
    filter(!is.na(!!response)) %>%
    mutate(diag_class = case_when(!!response == 1 & !!predictor == 1 ~ "TP",
                                  !!response == 1 & !!predictor == 0 ~ "FN",
                                  !!response == 0 & !!predictor == 0 ~ "TN",
                                  !!response == 0 & !!predictor == 1 ~ "FP")) %>% 
    select(diag_class, n) %>% 
    pivot_wider(names_from = "diag_class", values_from = "n", values_fill = list(n= 0)) %>% 
    mutate(Sensitivity = TP/(TP + FN),
           Specificity = TN/(TN + FP),
           PPV = TP/(TP+FP),
           NPV = TN/(FN+TN))
  
}

# Function to convert a set of dummy variables to a factor
# Args: data frame to reshape
# dummy_list = character vector or vars() list of the dummy variables to be converted
# id_var = variable identifying each row
# out_var = string giving name of output factor
# baseline_group = string giving factor level to be used as the baseline
DummyToFactor <- function(dat, dummy_list, id_var, out_var, baseline_group = NULL){
  # Select dataset to reshape
  selected_data <- if("quosures" %in% class(dummy_list)){
    {{dat}} %>% 
      select({{id_var}}, !!!dummy_list)
  } else {
    {{dat}} %>% 
      select({{id_var}}, !!dummy_list)
  }
  
  # Put into long format
  long_data <- {{dat}} %>% 
    select({{id_var}}) %>% 
    left_join(selected_data %>% 
                pivot_longer(cols = -{{id_var}}) %>% 
                filter(value == 1|value == T))
  
  # Set to NA ids with more than one dummy value recorded
  out_var_name <- enquo(out_var)
  long_data_filtered <- long_data %>% 
    count({{id_var}}) %>% 
    inner_join(long_data) %>% 
    mutate(name = if_else(n == 1, name, as.character(NA))) %>% 
    distinct() %>% 
    mutate_at("name", as.factor)
  
  long_data_formatted <- 
    if(!is.null(baseline_group)){
      long_data_filtered %>% 
        mutate(name = fct_relevel(name, baseline_group)) 
    } else {
      long_data_filtered
    } 
  
  long_data_formatted %>% 
    transmute({{id_var}}, !!out_var := name)
  
}
# DummyToFactor(assembled, vars(Smoke_Current, Smoke_Former, Smoke_Never), NISA_ID, "Smoke") 
# DummyToFactor(assembled, vars(Smoke_Current, Smoke_Former, Smoke_Never), NISA_ID, "Smoke", "Smoke_Never") 


# Function to calculate spherical equivalent
# Args: sphere = sphere in diopters
# cylinder = cylinder in diopters
CalculateSphericalEquivalent <- function(sphere, cylinder){
  if_else(cylinder == 0 | is.na(cylinder), sphere, sphere + cylinder/2)
}





