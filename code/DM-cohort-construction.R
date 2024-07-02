# Assemble cohort of all those with diabetes in the EMR dataset

library(tidyverse)
library(lubridate)

source(rprojroot::find_rstudio_root_file("code/EMR-helper-functions.R"))
source(rprojroot::find_rstudio_root_file("code/BIRAX-assemble-Medisoft-data.R"))
source(rprojroot::find_rstudio_root_file("code/BIRAX-assemble-NOA-data.R"))


# Find DMO eyes
DMO_eyes <- bind_rows(
  
  diagnoses %>% 
    filter(str_detect(DiagnosisDescription, "diabetic macular oedema"),
           !str_detect(DiagnosisDescription, "absent")) %>% 
    select(PatientID, EyeCode, DiagnosisDescription_DMO = DiagnosisDescription, DateofDiagnosis_DMO = DateofDiagnosis),
  
  clinical %>% 
    filter(str_detect(ClinicalTermDescription, "diabetic macular oedema"), 
           !str_detect(ClinicalTermDescription, "no diabetic macular oedema")) %>% 
    transmute(PatientID, EyeCode, DiagnosisDescription_DM = ClinicalTermDescription, DateofDiagnosis_DM = EncounterDate),
  
  dr_assessment %>% 
    select(PatientID, EyeCode, DiagnosisDescription_DMO = DRAssessmentDesc, DateofDiagnosis_DMO = EncounterDate) %>% 
    filter(str_detect(DiagnosisDescription_DMO, "DMO|oedema|CSMO"))
  
) %>% 
  group_by(PatientID, EyeCode) %>% 
  slice_min(order_by = DateofDiagnosis_DMO, with_ties = FALSE) %>% 
  ungroup()


# Find DR eyes
DR_eyes <- bind_rows(
  
  copath %>% 
    filter(CoPathologyCode %in% c("R317")) %>% 
    select(PatientID, EyeCode, DiagnosisDescription_DR = CoPathologyDesc, DateofDiagnosis_DR = EncounterDate),
  
  diagnoses %>% 
    filter(str_detect(DiagnosisDescription, "diabetic retinopathy"),
           !str_detect(DiagnosisDescription, "no diabetic retinopathy")) %>% 
    select(PatientID, EyeCode, DiagnosisDescription_DR = DiagnosisDescription, DateofDiagnosis_DR = DateofDiagnosis),
  
  clinical %>% 
    filter(str_detect(ClinicalTermDescription, "diabetic retinopathy"), 
           !str_detect(ClinicalTermDescription, "no diabetic retinopathy")) %>% 
    transmute(PatientID, EyeCode, DiagnosisDescription_DM = ClinicalTermDescription, DateofDiagnosis_DM = EncounterDate),
  
  dr_assessment %>% 
    select(PatientID, EyeCode, DiagnosisDescription_DR = DRAssessmentDesc, DateofDiagnosis_DR = EncounterDate) %>% 
    filter(str_detect(DiagnosisDescription_DR, "DMO|oedema|CSMO|none", negate = TRUE))
  
) %>% 
  group_by(PatientID, EyeCode) %>% 
  slice_min(order_by = DateofDiagnosis_DR, with_ties = FALSE) %>% 
  ungroup()





## Find full list of patients with either a systemic diabetes label, DMO or DR ##

# Full list of patients with DM and/or DR and/or DMO
# Note that there are some patients with DR or DMO but no systemic DM label
DM_patients <- 
  # Systemic DM label  
  diabetes %>% 
  filter(Diabetes) %>% 
  group_by(PatientID) %>% 
  slice_min(order_by = AgeDiagnosed, with_ties = FALSE) %>%
  ungroup() %>% 
  full_join(
    # DR status
    DR_eyes %>% 
      mutate(DR = 1) %>% 
      pivot_wider(id_cols =PatientID, names_from = EyeCode, values_from = DR, names_prefix = "DR_", values_fill = 0)
    , by = "PatientID") %>% 
  full_join(
    # DMO status
    DMO_eyes %>% 
      mutate(DMO = 1) %>% 
      pivot_wider(id_cols =PatientID, names_from = EyeCode, values_from = DMO, names_prefix = "DMO_", values_fill = 0) 
    , by = "PatientID") %>% 
  mutate(across(matches("DR|DMO"), ~if_else(is.na(.), 0, .)),
         DM = 1)


# Find first and last encounters for each patient in medical retina clinic
# Encounters are recorded at the patient level not the eye level
followup_dm <-
  encounters_raw[, .(
    FirstEncounterDate = EncounterDate[which.min(EncounterDateTime)],
    LastEncounterDate = EncounterDate[which.max(EncounterDateTime)],
    encounters = .N
  ), by = PatientID][DM_patients, on = "PatientID"][
    ,
    #  Calculate followup times
    # Note that some diagnoses were recorded prior to the first encounter.
    c("years_followup") :=
      list(
        interval(FirstEncounterDate, LastEncounterDate) / 
          dyears()
      )]



# For each visual acuity measurement, calculate the time since the first encounter (index_date)
va_raw_dm <- followup_dm[, .(PatientID, FirstEncounterDate)][ 
  # Join DM followup summary and visual_acuity tables
  visual_acuity, .(PatientID, EyeCode, FirstEncounterDate, EncounterDate, va_logmar, va_etdrs, va_category_snellen), 
  on = .(PatientID), nomatch = 0][
    # Calculate months since index date
    , months_since_index:=as.interval(FirstEncounterDate, EncounterDate)/dmonths()][
      # Calculate years since index date
      , years_since_index := as.interval(FirstEncounterDate, EncounterDate)/dyears()][  
        # Mark baseline measurements (closest measurement to baseline)  
        , baseline := abs(months_since_index) == min(abs(months_since_index)), by = .(PatientID)]



# List of candidate eyes
eye_raw_dm  <- 
  
  # All patients
  patients %>% 
  select(SiteID, PatientID, PerturbedDateofBirth, PerturbedCurrentAge, PerturbedAgeatDeath, Gender) %>% 
  inner_join(
    
    # All eyes potentially in the dataset (those not recorded in any table assume normal)
    expand_grid(PatientID = patients$PatientID, EyeCode = c("R", "L")), by = "PatientID") %>% 

  # DM patients
  left_join(DM_patients %>% 
              select(PatientID, AgeDiagnosed, DiabetesTypeDesc, DM), by = "PatientID") %>% 
  mutate(DM = if_else(is.na(DM), 0, DM)) %>% 
  
  # DR status of eyes
  left_join(DR_eyes %>% 
  select(PatientID, EyeCode, DateofDiagnosis_DR)
  , by = c("PatientID", "EyeCode")) %>% 
  mutate(DR = if_else(is.na(DateofDiagnosis_DR), 0, 1)) %>% 

  # DM status of eyes
  left_join(DMO_eyes %>% 
            select(PatientID, EyeCode, DateofDiagnosis_DMO)
          , by = c("PatientID", "EyeCode")) %>% 
  mutate(DMO = if_else(is.na(DateofDiagnosis_DMO), 0, 1)) %>% 

  # Calculate age at first encounter (index date)
  left_join(followup_dm %>% 
              select(PatientID, FirstEncounterDate, encounters, years_followup), by = c("PatientID")) %>% 
  mutate(index_age =  interval(PerturbedDateofBirth, FirstEncounterDate)/dyears()) %>% 
  mutate(years_observed = coalesce(PerturbedCurrentAge, PerturbedAgeatDeath) - index_age) %>% 
  
  # Join VA measured on the index date or up to 2 weeks prior as the baseline
  left_join(va_raw_dm %>%
              filter(baseline, months_since_index > -0.5 & months_since_index <=0) %>% 
              transmute(PatientID, EyeCode, va_logmar, "ETDRS letters" = va_etdrs),
            by = c("PatientID", "EyeCode")) %>% 
  
  # Count number of OCT thickness measurements as indicator of whether OCT took place
  left_join(oct_thickness[, .(oct_thickness_n = .N), by = .(PatientID, EyeCode)],
            by = c("PatientID", "EyeCode")) %>% 
  
  # Exclude eyes with no DM diagnosis at any point
  mutate(exclude_no_dm = DM == 0,
         
         
         # Those without a baseline VA measurement
         exclude_no_va = is.na(va_logmar),
         
         # Those without any OCT thickness measurements
         exclude_no_thickness = is.na(oct_thickness_n)
         
         
  )


eye_raw_dm %>% 
  count(across(matches("exclude")))

# Apply the exclusion critera
eye_dm <-
  eye_raw_dm %>% 
  filter(!exclude_no_dm,
         !exclude_no_va,
         !exclude_no_thickness)



# Just the VA history relating to the selected eyes
# Both before and after GA diagnosis
va_history_dm <- va_raw_dm %>% 
  inner_join(eye_dm %>% 
               select(PatientID, EyeCode), by = c("PatientID", "EyeCode")) %>% 
  # Roll forward baseline VA measurements so that all series start at the index date 
  # (max roll forward is 2 weeks)
  mutate(months_since_index = if_else(baseline & months_since_index > -0.5 & months_since_index < 0, 0, months_since_index))


