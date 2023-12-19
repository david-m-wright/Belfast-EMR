# Helper functions
library(janitor)
library(tidyverse)
library(skimr)

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



# Function to format intervals generated using the cut() function
# Args:
# cut_intervals: factor generated using cut()
# lab: character string to put between the two boundaries
# FormatIntervals <- function(cut_intervals, lab){
#   fct_relabel(cut_intervals, ~str_replace(., "([0-9]+)\\]", "<=\\1")) %>% 
#     fct_relabel(~str_replace(., "([0-9]+)\\)", "<\\1")) %>% 
#     fct_relabel(~str_replace(., "\\[([0-9]+)", "\\1<=")) %>% 
#     fct_relabel(~str_replace(., "\\(([0-9]+)", "\\1<")) %>% 
#     fct_relabel(~str_replace(., ",", lab))
# }



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

# Function to calculate mean and range for continuous variables
# Args: x = numeric vector to summarise
# output_n = logical indicating whether to return the number of non-missing cases in square brackets
# sig_digits = number of significant digits to report
SummariseMeanRange <- function(x, output_n = F, sig_digits = 1){
  mean_range <- paste0(signif(mean(x, na.rm = T), digits = sig_digits), 
                         " (", 
                         signif(min(x, na.rm = T), digits = sig_digits), 
                         ", ",
                         signif(max(x, na.rm = T), digits = sig_digits), 
                         ")")
  if(output_n) {
    paste0(mean_range, " [", 
           sum(!is.na(x)),
           "]")
  } else {
    mean_range
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
  
  # Note if the column variable is set to NULL
  enquo_var <- enquo(col_var)
  
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
                                                  `Mean (Range)` = SummariseMeanRange(.data[[.x]], sig_digits = sig_digits, output_n = output_n),
                                                  ` ` = paste0(sum(.data[[.x]]), " (100%)"),
                                                  SummariseSD(.data[[.x]], sig_digits = sig_digits, output_n = output_n)),
                                       .groups = "drop") %>%
                             mutate(Value = "Value") %>%
                             pivot_wider(names_from =   -n,  values_from = n) %>%
                             rename_with(~str_replace(., "_Value", ""), ends_with("_Value")) %>% 
                             mutate(Level = .y)
                         } else { 
                           
                           # For one way table
                           if(rlang::quo_is_null(enquo_var)){
                             out_tab <- count(data, .data[[.x]])
                           } else {
                             out_tab <- count(data, .data[[.x]], {{col_var}}) %>% 
                               pivot_wider(names_from = {{col_var}},  values_from = n, values_fill = 0, names_sort = T) 
                           }
                           # For two way table
                           out_tab %>% 
                             mutate(across(-1, ~paste0(., " (", signif(./sum(.)*100, digits = sig_digits), "%)"))) %>%
                             # rename_with(~"Value", .cols = (matches("^n$"))) %>%
                             # rename(Level = .data[[.x]]) %>%
                             rename(Level = all_of(.x)) %>%
                             mutate(across(Level, as.character))
                           
                           # DiscreteFrequency(data, row_var = .data[[.x]], col_var = {{col_var}}, sig_digits = sig_digits)
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


# Perform association tests for each of a set of continuaous variables
# Classified by a single discrete variable
# For numeric variables a two sample t-test will be performed
# For factors a chi-squared test will be performed
# var1_list = character vector or vars() specification listing variables to test (rows)
# var2 = unquoted name of variable to test by (columns)
SummariseAssocTest <-  function(dat, var1_list, var2){
  
  cross_tabs <- if("quosures" %in% class(var1_list)){
    {{var1_list}}
  } else {
    var1_list %>% 
      syms() 
  }
  
    dat %>% 
    select(all_of(var1_list)) %>% 
     map(~if(is.numeric(.)){
       t.test(. ~ pull(dat, {{var2}}))
     } else if(is.factor(.)) {
       stats::chisq.test(x = ., y = pull(dat, {{var2}}))
     }
     ) %>%
      map_dfr(broom::glance, .id = "Variable") %>%
                      transmute(Variable,
                                `Test stat` = round(statistic, digits = 2),
                                df = round(parameter, digits = 1),
                                P = format.pval(p.value, eps = 0.001, digits = 2))
}

 # SummariseAssocTest(oct_visits, noa_fluid, injected) 
# SummariseAssocTest(oct_visits, noa_grid, "injected")



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


# Function to display the diagnostic performance stats of a continuous predictor variable in AUC terms
# Args:
# dat = tibble containing data
# response = column containing response variable
# predictor_list = character vector or vars() specification listing variables to test (rows)
# Value: tibble containing the diagnostic stats
DiagnosticROC <- function(dat, response, predictor_list){
  # response = enquo(response)
  # predictor = enquo(predictor)
  cross_tabs <- if("quosures" %in% class(predictor_list)){
    {{predictor_list}}
  } else {
    predictor_list %>% 
      syms() 
  }
  
  dat %>% 
    select(all_of(predictor_list)) %>% 
    select(where(is.numeric)) %>% 
map(~ roc(response = pull(dat, {{response}}),  predictor = .x, ci = TRUE))  %>% 
  enframe(name = "Variable", value = "ROC") %>% 
  # AUCs 
  mutate(auc = map(.x = ROC, .f = ~ .x$ci),
         AUC = map_dbl(.x = ROC, .f = ~ .x$auc),
         lcl = map_dbl(.x = ROC, .f = ~.x$ci[1]),
         ucl = map_dbl(.x = ROC, .f = ~.x$ci[3]),
         across(c(AUC, lcl, ucl), formatC, format = "f", digits = 3),
         AUC_ci = paste0("(", lcl, ", ", ucl, ")")) %>% 
    select(Variable, AUC, AUC_ci)

}
# DiagnosticROC(oct_visits, injected, noa_fluid)

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


# Function to produce a circular plot of the ETDRS grid with labels corresponding to values in the source data 
# Args:
# dat = data.frame with joining columns "region" and "circ" 
# label_col = bare variable name containing the labels to apply to the grid
PlotETDRSGrid <- function(dat, label_col){
  
  label_value <- enquo(label_col)
  
  etdrs_rect <- tibble(region = as_factor(c("CS", "SUP", "NAS", "INF", "TEM", "SUP", "NAS", "INF", "TEM")),
                       prefix = c("", rep("i", 4), rep("o", 4)),
                       circ = c(1, rep(3, 4), rep(6, 4)),
                       xmin = c(0, rep(0.5, 4), rep(1.5, 4)),
                       xmax = c(0.5, rep(1.5, 4), rep(3, 4)),
                       ymin = c(0, 0, 0.25, 0.5, 0.75, 0, 0.25, 0.5, 0.75),
                       ymax = c(1, 0.25, 0.5, 0.75, 1, 0.25, 0.5, 0.75,1),
                       xtext = c(0, xmax[-1] - (xmax[-1]-xmin[-1])/2),
                       ytext = c(0, ymax[-1] - (ymax[-1]-ymin[-1])/2))
  
  dat %>% 
    inner_join(etdrs_rect, by = c("region", "circ")) %>% 
    
    ggplot(aes(x = xtext, y = ytext)) +
    geom_rect(data = filter(etdrs_rect, region != "CS"), aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "white", colour = "black", linewidth = 1.5) + 
    geom_rect(data = filter(etdrs_rect, region == "CS"), aes(xmin = xmin+0.01, xmax = xmax-0.01, ymin = ymin+0.01, ymax = ymax-0.01), fill = "white") +
    geom_text(aes(label = !!label_value)) +
    scale_y_continuous(breaks = NULL) +
    scale_x_continuous(breaks = NULL, sec.axis = sec_axis(trans = ~., breaks = NULL, name="Nasal")) +
    theme_light() +
    theme(panel.border = element_blank(), panel.grid = element_blank()) +
    labs(x= "Temporal", y = NULL)+
    coord_polar(theta = "y", start = 315*pi/180)
  
}
# 
# baseline_fluid_grid %>% 
#   filter(str_detect(abbr_metric, "IRF Vol")) %>% 
#   PlotETDRSGrid(valform) 



# Custom skimmers for dealing with numeric values
# Mean and sd
skim_mean_sd <- function(x, sig_digits = 1){
  paste0(signif(mean(x, na.rm = TRUE), digits = sig_digits), 
         " (", 
         signif(sd(x, na.rm = TRUE), digits = sig_digits),
         ")")
}
skim_mean_sd(1:10)

# Median and range
skim_median_range <- function(x, sig_digits = 1){
  paste0(signif(median(x, na.rm = T), digits = sig_digits), 
         " (", 
         signif(min(x, na.rm = T), digits = sig_digits), 
         ", ",
         signif(max(x, na.rm = T), digits = sig_digits), 
         ")")
}
skim_median_range(1:10) 

# Median and IQR
skim_median_iqr <- function(x, sig_digits = 1){
  paste0(signif(median(x, na.rm = T), digits = sig_digits), 
         " (", 
         signif(quantile(x, probs = c(0.25), na.rm = T), digits = sig_digits), 
         ", ",
         signif(quantile(x, probs = c(0.75), na.rm = T), digits = sig_digits), 
         ")")
}
skim_median_iqr(1:10) 

# Custom skimmer list for numeric values
skim_numeric <- skim_with(base = list(n = ~ as.character(n_complete(.))),
                          numeric = list(
                            `Mean (SD)` = ~ skim_mean_sd(., sig_digits = 3),
                            `Median (Range)` = ~skim_median_range(., sig_digits = 3),
                            `Median (IQR)` = ~skim_median_iqr(., sig_digits = 3),
                            mean = NULL,
                            sd = NULL,
                            p0  = NULL,
                            p25 = NULL,
                            p50  = NULL,
                            p75 = NULL,
                            p100 = NULL,
                            hist = NULL)
)

# Function to format the output from skim_numeric()
format_skim_numeric <- function(x, col_var) {x %>% 
    yank("numeric") %>% 
    rename(Variable = "skim_variable") %>% 
    pivot_longer(cols = c(n, `Mean (SD)`, `Median (Range)`, `Median (IQR)`), names_to = "Level") %>% 
    pivot_wider(names_from = {{col_var}}, values_from = value) 
}
