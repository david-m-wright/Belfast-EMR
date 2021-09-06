# Helper functions
library(janitor)
library(tidyverse)


# Function to reformat the intervals produced by cut() into more familiar inequality notation
# Args: x = character vector to reformat,
# unit = character string to be placed in the middle of the inequality
IntervalToInequality <- function(x, unit = "x"){
  enframe(x) %>% 
    separate(value, sep = ",", into = c("lower", "upper"), remove = FALSE) %>% 
    mutate(
      # lclosed = if_else(str_detect(lower, "\\["), "\u2264", "<"),
      # uclosed = if_else(str_detect(upper, "\\]"), "\u2264", "<"),
      lclosed = if_else(str_detect(lower, "\\["), "<=", "<"),
      uclosed = if_else(str_detect(upper, "\\]"), "<=", "<"),
      
      across(c(lower, upper), ~str_remove(., "(\\[|\\(|\\)|\\])")),
      output = if_else(!is.na(value), paste(lower, lclosed, unit, uclosed, upper), as.character(NA))) %>% 
    pull(output)
}
# IntervalToInequality(baseline$Myopia)


# Function to select most frequently occuring non-missing value from a set of values
# Or first (using sort()) in case of ties
# Args:
# x = vector of values to count votes in
# Value: value with the most votes
MaxVotes <- function(x){
  if(all(is.na(x))){return(NA)}
  names(sort(table(x), decreasing = T))[1]
}


# Function to calculate SD and mean for continous variables
# Args: x = numeric vector to summarise
# output_n = logical indicating whether to return the number of non-missing cases in square brackets
# sig_digits = number of significant digits to report
SummariseSD <- function(x, output_n = F, sig_digits = 1){
  mean_sd <- paste0(signif(mean(x, na.rm = T), digits = sig_digits), 
                    " (", 
                    signif(sd(x, na.rm = T), digits = sig_digits), 
                    ")")
  if(output_n){
    paste0(mean_sd, " [", 
           sum(!is.na(x)),
           "]") 
  } else {
    mean_sd
  }
}

# Function to calculate median and interquartile range for continuous variables
# Args: x = numeric vector to summarise
# output_n = logical indicating whether to return the number of non-missing cases in square brackets
# sig_digits = number of digits to report
SummariseMedianIQR <- function(x, output_n = F, sig_digits = 1){
  median_iqr <- paste0(signif(median(x, na.rm = T), digits = sig_digits), 
                       " (", 
                       signif(quantile(x, probs = c(0.25), na.rm = T), digits = sig_digits), 
                       ", ",
                       signif(quantile(x, probs = c(0.75), na.rm = T), digits = sig_digits), 
                       ")")
  if(output_n) {
    paste0(median_iqr, " [", 
           sum(!is.na(x)),
           "]")
  } else {
    median_iqr
  }
}


# Function to calculate median and range for continuous variables
# Args: x = numeric vector to summarise
# output_n = logical indicating whether to return the number of non-missing cases in square brackets
# sig_digits = number of significant digits to report
SummariseMedianRange <- function(x, output_n = F, sig_digits = 1){
  median_range <- paste0(signif(median(x, na.rm = T), digits = sig_digits), 
                         " (", 
                         signif(min(x, na.rm = T), digits = sig_digits), 
                         ", ",
                         signif(max(x, na.rm = T), digits = sig_digits), 
                         ")")
  if(output_n) {
    paste0(median_range, " [", 
           sum(!is.na(x)),
           "]")
  } else {
    median_range
  }
}

# Function to generate a frequency table for a discrete variable (n, %),
# optionally cross tabulated by a second variable
# Args: data = data.frame containing data to tabulate
# row_var = unquoted name of variable to tabulate
# col_var = for a two-way table, unquoted name of variable to cross tabulate
# sig_digits = integer giving number of significant digits to report
# Value: tibble containing formatted counts and column percentages
DiscreteFrequency <- function(data, row_var, col_var, sig_digits = 3){
  count(data, {{row_var}}, {{col_var}}) %>% 
    pivot_wider(names_from = {{col_var}},  values_from = n, values_fill = 0, names_sort = T) %>%
    mutate(across(-1, ~paste0(., " (", signif(./sum(.)*100, digits = sig_digits), "%)"))) %>% 
    rename_with(~"Value", .cols = (matches("^n$"))) %>% 
    rename(Level = {{row_var}}) %>%
    mutate(across(Level, as.character))
}

#  
#  study_eye_data %>% 
#    DiscreteFrequency(row_var = DM_Type)
# # 
#  study_eye_data %>% 
#    DiscreteFrequency(DM_Type, col_var = DM_DR_status)
#  


# Function to generate a table of descriptive statistics,
# displaying summaries of continuous variables (e.g. mean, SD) or distribution of discrete variables (e.g. %) as appropriate
# Args: dat = data.frame containing data to tabulate
# Each column will form a row of the descriptive table
# col_var = for a two-way table, unquoted name of variable to cross tabulate
# type = named character vector specifying the type of summary statistic for continuous variables
# default is mean and SD, can also specify median and IQR or median and range
# sig_digits = integer indicating the number of significant digits to report for each summary
# output_n = logical indicating whether to return the number of non-missing cases in square brackets
# Value: tibble containing formatted table of descriptives
GenerateDescriptives <- function(data, col_var = NULL, type = NULL, sig_digits = 3, output_n = FALSE){
  
  # Columns to tabulate
  tab_data <- data %>% 
    select(-{{col_var}})
  
  # Select type of numerical summary for continuous variables 
  ctype <- setNames(rep("Mean (SD)", length.out = ncol(tab_data)),
                    map(names(tab_data), rlang::as_name))
  ctype[match(names(type), names(ctype))] <- type
  
  # Cross tabulate each variable in the list
  base_table <- map2_dfr(.x = setNames(names(ctype), names(ctype)), 
                         .y = ctype, 
                         .f = ~if(is.numeric(data[[.x]]))
                         { 
                           data %>%
                             group_by({{col_var}}) %>%
                             summarise(n = switch(.y,
                                                  `Median (IQR)` = SummariseMedianIQR(.data[[.x]], sig_digits = sig_digits, output_n = output_n),
                                                  `Median (Range)` = SummariseMedianRange(.data[[.x]], sig_digits = sig_digits, output_n = output_n),
                                                  SummariseSD(.data[[.x]], sig_digits = sig_digits, output_n = output_n)),
                                       .groups = "drop") %>%
                             mutate(Value = "Value") %>%
                             pivot_wider(names_from =   -n,  values_from = n) %>%
                             rename_with(~str_replace(., "_Value", ""), ends_with("_Value")) %>% 
                             mutate(Level = .y)
                         } else { 
                           DiscreteFrequency(data, row_var = .data[[.x]], col_var = {{col_var}}, sig_digits = sig_digits)
                         }, .id = "Variable")
  
  base_table %>% 
    select(Variable, Level, everything())
}  

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



