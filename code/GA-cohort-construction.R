# Assemble GA and fibrosis project cohorts

library(tidyverse)
library(lubridate)

source(rprojroot::find_rstudio_root_file("code/EMR-helper-functions.R"))
source(rprojroot::find_rstudio_root_file("code/BIRAX-assemble-Medisoft-data.R"))
source(rprojroot::find_rstudio_root_file("code/BIRAX-assemble-NOA-data.R"))



# Find eyes with GA and find first GA diagnosis date for each eye (index date)
diagnosis_ga <- 
  bind_rows(
  # Copathology table
  copath %>% 
  filter(CoPathologyCode == "R338") %>% 
  select(PatientID, EyeCode, DiagnosisDescription_GA = CoPathologyDesc, DateofDiagnosis_GA = EncounterDate),
  # Diagnoses table  
  diagnoses %>% 
      filter(str_detect(DiagnosisDescription, "geographic atrophy"), 
             !str_detect(DiagnosisDescription, "non geographic atrophy")) %>% 
      select(PatientID, EyeCode, DiagnosisDescription_GA = DiagnosisDescription, DateofDiagnosis_GA = DateofDiagnosis),
  
  # Clincial findings table
  clinical %>% 
    filter(str_detect(ClinicalTermDescription, "geographic atrophy"), 
           !str_detect(ClinicalTermDescription, "non geographic atrophy")) %>% 
    transmute(PatientID, EyeCode, DiagnosisDescription_GA = ClinicalTermDescription, DateofDiagnosis_GA = EncounterDate)

  ) %>% 
  mutate(DateofDiagnosis_GA = as.POSIXct(as.Date(DateofDiagnosis_GA)))  %>% 

  group_by(PatientID, EyeCode) %>% 
  slice_min(order_by =  DateofDiagnosis_GA, with_ties = FALSE) %>% 
  ungroup() %>% 
  as.data.table()


# diagnosis_ga
# diagnosis_ga %>% 
#   distinct(PatientID)  

# Find first and last encounters for each patient in medical retina clinic
followup_ga <-
  encounters_raw[, .(
    FirstEncounterDate = EncounterDate[which.min(EncounterDateTime)],
    LastEncounterDate = EncounterDate[which.max(EncounterDateTime)],
    encounters = .N
  ), by = PatientID][diagnosis_ga, on = "PatientID"][
    ,
                                                #  Calculate followup times
                                                # Note that some diagnoses were recorded prior to the first encounter.
                                                c("years_followup",
                                                  "years_pre_GA",
                                                  "years_with_GA") :=
                                                  list(
                                                    interval(FirstEncounterDate, LastEncounterDate) / 
                                                      dyears(),
                                                    interval(FirstEncounterDate, DateofDiagnosis_GA) /
                                                      dyears(),
                                                    interval(DateofDiagnosis_GA, LastEncounterDate) /
                                                      dyears()
                                                  )]

# followup_ga
# summary(followup_ga$years_with_GA)
# followup_ga %>% filter(years_with_GA == 0)
# #  Small number diagnosed before first encounter
# followup_ga %>% filter(years_pre_GA<0)
# encounters_raw[PatientID == "065741D6-A4C7-A4C5-4176-297FC352993C"]


# For each visual acuity measurement, calculate the time since the GA diagnosis (index_date)
va_raw_ga <- followup_ga[, .(PatientID, EyeCode, DateofDiagnosis_GA)][ 
  # Join GA followup summary and visual_acuity tables
  visual_acuity, .(PatientID, EyeCode, DateofDiagnosis_GA, EncounterDate, va_logmar, va_etdrs, va_category_snellen), 
  on = .(PatientID, EyeCode), nomatch = 0][
    # Calculate months since index date
    , months_since_index:=as.interval(DateofDiagnosis_GA, EncounterDate)/dmonths()][
      # Calculate years since index date
      , years_since_index := as.interval(DateofDiagnosis_GA, EncounterDate)/dyears()][  
        # Mark baseline measurements (closest measurement to baseline)  
        , baseline := abs(months_since_index) == min(abs(months_since_index)), by = .(PatientID, EyeCode)]



eye_raw_ga  <- 

  # All patients
  patients %>% 
  select(SiteID, PatientID, PerturbedDateofBirth, PerturbedCurrentAge, PerturbedAgeatDeath, Gender) %>% 
  inner_join(
  
  # All eyes potentially in the dataset (those not recorded in any table assume normal)
  expand_grid(PatientID = patients$PatientID, EyeCode = c("R", "L")), by = "PatientID") %>% 
    
  #   # All eyes with a diagnosis (some eyes not recorded at all so assume normal)
  # diagnoses %>% 
  # distinct(PatientID, EyeCode)
  # , by = "PatientID") 
 

  # Add diagnoses for DM (see definition from DMO cohort script)
  
  # Calculate age at first GA diagnosis (index date)
  left_join(followup_ga, by = c("PatientID", "EyeCode")) %>% 
  mutate(index_age =  interval(PerturbedDateofBirth, DateofDiagnosis_GA)/dyears()) %>% 
  mutate(years_observed = coalesce(PerturbedCurrentAge, PerturbedAgeatDeath) - index_age) %>% 

  # Join VA measured on the index date or up to 2 weeks prior as the baseline
  left_join(va_raw_ga %>%
              filter(baseline, months_since_index > -0.5 & months_since_index <=0) %>% 
              transmute(PatientID, EyeCode, va_logmar, "ETDRS letters" = va_etdrs),
            by = c("PatientID", "EyeCode")) %>% 

  # Count number of OCT thickness measurements as indicator of whether OCT took place
  left_join(oct_thickness[, .(oct_thickness_n = .N), by = .(PatientID, EyeCode)],
            by = c("PatientID", "EyeCode")) %>% 

  # Exclude eyes with no GA diagnosis at any point
  mutate(exclude_no_GA = is.na(DateofDiagnosis_GA),
         
         
         # Those without a baseline VA measurement
         exclude_no_va = is.na(va_logmar),
         
         # Those without any OCT thickness measurements
         exclude_no_thickness = is.na(oct_thickness_n)
         
         
  )


eye_raw_ga %>% 
  count(across(matches("exclude")))
  
# Apply the exclusion critera
eye_ga <-
  eye_raw_ga %>% 
  filter(!exclude_no_GA,
         !exclude_no_va,
         !exclude_no_thickness)



# Just the VA history relating to the selected eyes
# Both before and after GA diagnosis
va_history_ga <- va_raw_ga %>% 
  inner_join(eye_ga %>% 
               select(PatientID, EyeCode), by = c("PatientID", "EyeCode")) %>% 
  # Roll forward baseline VA measurements so that all series start at the index date 
  # (max roll forward is 2 weeks)
  mutate(months_since_index = if_else(baseline & months_since_index > -0.5 & months_since_index < 0, 0, months_since_index))


