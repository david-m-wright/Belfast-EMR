# Descriptive statistics for BIRAX dataset

library(skimr)
library(tidyverse)
library(data.table)
source(rprojroot::find_rstudio_root_file("code/BIRAX-assemble-cohort.R"))


# Timing of functional response to treatment

# Classifying treatment response according to VA measurements
# Three stages to response
# Primary = first VA visit after third injection but must be <12 months after first injection
# Secondary = second VA visit after the third injection but must be <12 months after first injection
# Late = first VA visit >= 12 months after first injection

# Find date of third injection
# Add join_date column to prevent confusion about which date is reported after join
injections_third <- as.data.table(injections)[,join_date:=EncounterDate,][
    injection_number == 3, .(PatientID, EyeCode, injection_date = EncounterDate, join_date, injection_number)]

# Join the third injections to the VA history
va_response_raw <- injections_third[va_history[, join_date := EncounterDate],
  roll = Inf, on = .(PatientID, EyeCode, join_date)
][# Number the VA measurements that took place after 3rd injection
  , va_seq := row_number(EncounterDate), by = .(PatientID, EyeCode, injection_number)
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

# # Write out the VA response data frame to spss
# haven::write_sav(va_response, "//fas8200main-n2/OphBelfast/va_response.sav")
# 
# va_history %>% 
#   filter(snapshot,
#          follow_up_year != 0 , follow_up_year != 0.5) %>% 
#   haven::write_sav("//fas8200main-n2/OphBelfast/va_follow_up.sav")


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



## Fluid characteristics by response phase

# Assign each fluid measurement to a phase
# This uses the same phase event timings for each eye as the VA analysis


# Join the third injections to the fluid history
fluid_response_raw <- injections_third[fluid_history[, join_date := as.POSIXct(OCTDate)],
                                    roll = Inf, on = .(PatientID, EyeCode, join_date)
][# Number the fluid measurements that took place after 3rd injection
  , fluid_seq := row_number(OCTDate), by = .(PatientID, EyeCode, injection_number)
][
  # Indicate first fluid measurement after 3rd injection (only if <12 months after first injection)
  , primary := injection_number == 3 & fluid_seq == 1 & months_since_index < 12
][
  # Indicate second fluid measurement after 3rd injection (only if <12 months after first injection)
  , secondary := injection_number == 3 & fluid_seq == 2 & months_since_index < 12
][
  # Indicate first VA measurement >=12 months after first injection (in second treatment year)
  , late := injection_number == 3 & .I == .I[which.min(fluid_seq)] & treatment_year == 2
  , by = .(PatientID, EyeCode, treatment_year)
][
  # Combine the three classifications into a single variable
  , response := factor(case_when(primary ~ "Primary",
                                 secondary ~ "Secondary",
                                 late ~ "Late"), levels = c("Primary", "Secondary", "Late"))
] 
fluid_response_raw %>% count(response)
# Calculate proportional change since baseline
setkey(fluid_response_raw, PatientID, EyeCode, OCTDate)
fluid_response_raw[, TRFVolumeNl_change_prop := TRFVolumeNl/TRFVolumeNl[1], by = .(PatientID, EyeCode)]

# Classify fluid response (total fluid) as good partial or non-response
fluid_response_raw[,
  fl_response := factor(case_when(TRFVolumeNl == 0 ~ "Good",
                                  TRFVolumeNl > 0 & TRFVolumeNl_change_prop <= 0.7 ~ "Partial",
                                  TRFVolumeNl > 0 & TRFVolumeNl_change_prop > 0.7 ~ "Non-responder"
  ),
  levels= c("Good", "Partial", "Non-responder"))
]

# Retain a single fluid measurement for each phase for each eye (not all eyes have all phases)
fluid_response <- fluid_response_raw[!is.na(response), ]


# fluid_response %>% 
  fluid_response_raw %>% 
  # mutate(fl_response = factor(case_when(TRFVolumeNl == 0 ~ "Good",
  #                                TRFVolumeNl > 0 & TRFVolumeNl_change_prop <= 0.7 ~ "Partial",
  #                                TRFVolumeNl > 0 & TRFVolumeNl_change_prop > 0.7 ~ "Non-responder"
  #                                ),
  #                             levels= c("Good", "Partial", "Non-responder"))) %>% 
  count(response, fl_response)
         
fluid_response_raw %>% filter(is.na(response)) %>% 
  select(PatientID, EyeCode, fl_response, TRFVolumeNl, TRFVolumeNl_change_prop)
  
fluid_response_raw[PatientID == 3 & EyeCode == "OD",.(PatientID, EyeCode,OCTDate, months_since_index, TRFVolumeNl, TRFVolumeNl_change_prop)]
fluid_response[PatientID == 3 & EyeCode == "OD",.(PatientID, EyeCode,OCTDate, months_since_index, TRFVolumeNl, TRFVolumeNl_change_prop, fl_response)]


# Classifying fluid changes by phase of functional response
fluid_phases <- fluid_response %>% 
  select(PatientID, EyeCode, response, fl_response) %>% 
  pivot_wider(values_from = fl_response, names_from = response) %>% 
  filter(!is.na(Primary) | !is.na(Secondary) | !is.na(Late))


# Sub-analysis by initial fluid response

# Find the time and number of injections required to achieve fluid resolution
# Join schematic: injections[fluid_response_raw, roll = TRUE] 
# Find the sequence number of the injection immediately before the resolution (rolling join)
first_resolution <- as.data.table(injections)[, join_date := EncounterDate][
    # Find the first fluid measurement at which the fluid resolution was achieved
    fluid_response_raw[fl_response == "Good" & !baseline, .SD[which.min(months_since_index)], by = .(PatientID, EyeCode)][
      , join_date := as.POSIXct(OCTDate)],  
    # Perform the rolling join
    roll = TRUE, on = .(PatientID, EyeCode, join_date)]
  

  


# Analysis by number of annual injections

# Number of injection year 1
injections_yr1 <- injections %>% 
  group_by(PatientID, EyeCode) %>% 
  filter(treatment_year == 1) %>% 
  summarise(injections_yr1 = n(), .groups = "drop")

# Number of injections per year aggregated across years 2- 3
injections_per_yr <- injections %>% 
  group_by(PatientID, EyeCode) %>% 
  filter(treatment_year == 2 | treatment_year == 3) %>% 
  summarise(injections_per_yr = n()/2, .groups = "drop") %>% 
  mutate(injections_per_yr_grp = cut(injections_per_yr, breaks = c(min(injections_per_yr), 4, 7, max(injections_per_yr)), 
                                     include.lowest = TRUE),
         injections_per_yr_grp = fct_relabel(injections_per_yr_grp, IntervalToInequality))

