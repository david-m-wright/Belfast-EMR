# Descriptive statistics for TASMC dataset

library(skimr)
library(tidyverse)
library(data.table)
source(rprojroot::find_rstudio_root_file("code/TASMC-assemble-cohort.R"))


# Timing of functional response to treatment

# Classifying treatment response according to VA measurements
# Three stages to response
# Primary = first VA visit after third injection but must be <12 months after first injection
# Secondary = second VA visit after the third injection but must be <12 months after first injection
# Late = first VA visit >= 12 months after first injection

# Find date of third injection
# Add join_date column to prevent confusion about which date is reported after join
injections_third <- as.data.table(injections)[,join_date:=Date,][
    injection_number == 3, .(PatientID, EyeCode, injection_date = Date, join_date, injection_number)]

# Join the third injections to the VA history
va_response_raw <- injections_third[va_history[, join_date := Date],
  roll = Inf, on = .(PatientID, EyeCode, join_date)
][# Number the VA measurements that took place after 3rd injection
  , va_seq := row_number(Date), by = .(PatientID, EyeCode, injection_number)
  ][
  # Indicate first VA measurement after 3rd injection (only if <12 months after first injection)
  , primary := injection_number == 3 & va_seq == 1 & months_since_index < 12
][
  # Indicate second VA measurement after 3rd injection (only if <12 months after first injection)
  , secondary := injection_number == 3 & va_seq == 2 & months_since_index < 12
][
  # Indicate first VA measurement >=12 months after first injection (in second treatment year)
  , late := injection_number == 3 & .I == .I[which.min(va_seq)] & treatment_year == 2
  , by = .(PatientID, EyeCode, treatment_year)
][
  # Combine the three classifications into a single variable
  , response := factor(case_when(primary ~ "Primary",
                         secondary ~ "Secondary",
                         late ~ "Late"), levels = c("Primary", "Secondary", "Late"))
] 


# Retain a single VA measurement for each phase for each eye (not all eyes have all phases)
va_response <- va_response_raw[!is.na(response)]


va_response_raw[PatientID == 3 & EyeCode == "OS"]
va_response_raw %>% count(primary, secondary, late)
va_response_raw %>% count(PatientID, EyeCode, response) %>% count(response, n)


va_response[PatientID == 3]
va_response %>% 
  count(PatientID, EyeCode) %>% count(n)


# Classifying visual acuity changes by phase of functional response
va_phases <- va_response %>% 
  select(PatientID, EyeCode, response, va_change_cat) %>% 
  pivot_wider(values_from = va_change_cat, names_from = response) %>% 
  filter(!is.na(Primary) | !is.na(Secondary) | !is.na(Late))



