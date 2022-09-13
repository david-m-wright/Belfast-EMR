# Assemble BIRAX analysis cohort (eye level)

library(tidyverse)
library(lubridate)

source(rprojroot::find_rstudio_root_file("code/BIRAX-assemble-Medisoft-data.R"))
source(rprojroot::find_rstudio_root_file("code/BIRAX-assemble-NOA-data.R"))


# Find first AMD diagnosis date for each eye
AMD_patients <- copath %>% 
  filter(CoPathologyCode %in% c("R316", "R337")) %>% 
  select(PatientID, EyeCode, DiagnosisDescription_AMD = CoPathologyDesc, DateofDiagnosis_AMD = EncounterDate) %>% 
  bind_rows(
    diagnoses %>% 
      filter(str_detect(DiagnosisDescription, "AMD|age-related macular degeneration"), 
             !str_detect(DiagnosisDescription, "suspected|no evidence|dry"),
             !is.na(DiagnosisICD10Code)) %>% 
      select(PatientID, EyeCode, DiagnosisDescription_AMD = DiagnosisDescription, DateofDiagnosis_AMD = DateofDiagnosis)) %>% 
  group_by(PatientID, EyeCode) %>% 
  slice_min(order_by =  DateofDiagnosis_AMD, with_ties = FALSE) %>% 
  ungroup()
# Note small number with dry AMD at first diagnosis.
# AMD_patients %>% count(DiagnosisDescription_AMD)

# Find RVO patients
RVO_patients <- copath %>% 
  filter(CoPathologyCode == "R318") %>% 
  select(PatientID, EyeCode, DiagnosisDescription_RVO = CoPathologyDesc, DateofDiagnosis_RVO = EncounterDate) %>% 
  bind_rows(
    diagnoses %>% 
      filter(str_detect(DiagnosisDescription, "occlusion")) %>% 
      select(PatientID, EyeCode, DiagnosisDescription_RVO = DiagnosisDescription, DateofDiagnosis_RVO = DateofDiagnosis)) %>% 
  group_by(PatientID, EyeCode) %>% 
  slice_min(order_by = DateofDiagnosis_RVO, with_ties = FALSE) %>% 
  ungroup()

# Find DMO/DR patients
DMO_DR_patients <- bind_rows(
  copath %>% 
    filter(CoPathologyCode %in% c("R318")) %>% 
    select(PatientID, EyeCode, DiagnosisDescription_DR = CoPathologyDesc, DateofDiagnosis_DR = EncounterDate),
  
  diagnoses %>% 
    filter(str_detect(DiagnosisDescription, "diabetic macular oedema"),
           !str_detect(DiagnosisDescription, "absent")) %>% 
    select(PatientID, EyeCode, DiagnosisDescription_DR = DiagnosisDescription, DateofDiagnosis_DR = DateofDiagnosis),
  
  dr_assessment %>% 
    select(PatientID, EyeCode, DiagnosisDescription_DR = DRAssessmentDesc, DateofDiagnosis_DR = EncounterDate) %>% 
    filter(DiagnosisDescription_DR != "none")
  
) %>% 
  group_by(PatientID, EyeCode) %>% 
  slice_min(order_by = DateofDiagnosis_DR, with_ties = FALSE) %>% 
  ungroup()



### Prepare series of measurements over time ###

# These analyses are all based on the injection history for each eye, starting at the date of the first injection

# Injection history summary by eye
# Find date of first injection (IndexDate) and final injection
# and number of injections for each eye
injection_summary_eye  <- injections_clean %>% 
  group_by(PatientID, EyeCode) %>% 
  summarise(index_date = min(EncounterDate),
            final_injection_date = max(EncounterDate), 
            total_injections = n(),
            .groups = "drop") %>% 
  mutate(total_intervals = total_injections - 1)


# For each visual acuity measurement, calculate the time since the first injection (index_date)
va_raw <- as.data.table(injection_summary_eye)[, .(PatientID, EyeCode, index_date)][ 
  # Join injection_summary_eye and visual_acuity tables
  as.data.table(visual_acuity), .(PatientID, EyeCode, index_date, EncounterDate, va_logmar, va_etdrs, va_category_snellen), on = .(PatientID, EyeCode)][
  # Calculate months since index date
  , months_since_index:=as.interval(index_date, EncounterDate)/dmonths()][
  # Mark baseline measurements (closest measurement to baseline)  
  , baseline := abs(months_since_index) == min(abs(months_since_index)), by = .(PatientID, EyeCode)]


## OCT thickness metrics
oct_thickness_raw <- as.data.table(injection_summary_eye)[
  , .(PatientID, EyeCode, index_date)][
  # Join injection_summary_eye and oct_thickness tables
  as.data.table(oct_thickness), on = .(PatientID, EyeCode)][
  # Calculate months since index date
  , months_since_index:=as.interval(index_date, ExamDate)/dmonths()][
  # Mark baseline measurements (closest measurement to baseline)  
  , baseline := abs(months_since_index) == min(abs(months_since_index)), by = .(PatientID, EyeCode)]

# Calculate changes in thickness measurements from baseline 
# Could do this with a set of baseline columns and a set of change columns
# ss <- oct_thickness_raw %>% 
#   .[, .(PatientID, EyeCode, ExamDate, OuterSuperior_Volume, OuterInferior_Volume)] 
#   cnames <- c("OuterSuperior_Volume", "OuterInferior_Volume")
# Change x[1] to x[baseline]
#   ss[, .SDcols = cnames,  (paste0(cnames, "_change")) := lapply(.SD, function(x) x - x[1]), by=list(PatientID, EyeCode)] %>% 
#     print()


## NOA fluid measurements
  # For each OCT thickness measurement, calculate the time since the first injection (index_date)
fluid_raw <- as.data.table(injection_summary_eye)[
    , .(PatientID, EyeCode, index_date)][
    # Join injection_summary_eye and oct_thickness tables
    as.data.table(noa), on = .(PatientID, EyeCode)][
    # Calculate months since index date
    , months_since_index:=as.interval(index_date, OCTDate)/dmonths()][
    # Mark baseline measurements (closest measurement to baseline)  
    , baseline := abs(months_since_index) == min(abs(months_since_index)), by = .(PatientID, EyeCode)]


### Assemble eye level dataset (single row per eye) ###

eye_raw <- patients %>% 
  
  # Add diagnoses for AMD, RVO and DMO/DR
  left_join(AMD_patients, by = "PatientID") %>% 
  left_join(RVO_patients, by = c("PatientID", "EyeCode")) %>% 
  left_join(DMO_DR_patients, by = c("PatientID", "EyeCode")) %>% 
  
  # Calculate age at first injection (index date)
  left_join(injection_summary_eye, by = c("PatientID", "EyeCode")) %>% 
  mutate(index_age =  as.interval(PerturbedDateofBirth, index_date)/dyears()) %>% 
  # Calculate years observed (index to current age/death)
  mutate(years_observed = coalesce(PerturbedCurrentAge, PerturbedAgeatDeath) - index_age) %>% 
  # Calculate years treated (index to final injection)
  mutate(years_treated = interval(index_date, final_injection_date)/dyears()) %>% 
  
  # Join VA measured on the index date or up to 2 weeks prior as the baseline
  left_join(va_raw %>%
              filter(baseline, months_since_index > -0.5 & months_since_index <=0) %>% 
              transmute(PatientID, EyeCode, EncounterDate, va_logmar, "ETDRS letters" = va_etdrs),
            by = c("PatientID", "EyeCode")) %>% 
  
  # Join presence of OCT thickness measurements (closest to index date but up to 2 weeks prior)
  # Even though an eye may have no baseline it may have later measurements that can be included
  # although change metrics cannot be calculated for these.
  left_join(oct_thickness_raw %>% 
              filter(baseline, months_since_index > -0.5 & months_since_index <=0) %>% 
              transmute(PatientID, EyeCode, thickness_measurement = baseline), 
            by = c("PatientID", "EyeCode")) %>% 
  
  # Join presence of fluid measurements (closest to index date but upto 2 weeks prior)
  left_join(fluid_raw %>% 
              filter(baseline, months_since_index > -0.5 & months_since_index <=0) %>% 
              transmute(PatientID, EyeCode, fluid_measurement = baseline), 
            by = c("PatientID", "EyeCode")) %>% 
  
  mutate(across(c(thickness_measurement, fluid_measurement), ~if_else(is.na(.), FALSE, .))) %>%
  
  # Exclusion criteria
  mutate(
    # 50 or older on date of first injection (index date)
    exclude_age = index_age < 50 |is.na(index_age),
    
    # Exclude no AMD diagnosis
    exclude_no_AMD = is.na(DiagnosisDescription_AMD),
    
    # Exclude <3 anti-VEGF injections
    exclude_lt3_injections = if_else(total_injections < 3 | is.na(total_injections), TRUE, FALSE),
    
    # Diagnosis of DMO/DR or RVO
    exclude_DR_DMO = !is.na(DiagnosisDescription_DR),
    exclude_RVO = !is.na(DiagnosisDescription_RVO),
    
    # No injections after index date (so no treatment intervals can be calculated)
    exclude_no_intervals = total_intervals == 0,
    
    # Those without a baseline VA measurement
    exclude_no_va = is.na(va_logmar)) 


# Apply the exclusion criteria
eye <- eye_raw %>% 
  filter(!exclude_no_AMD, !exclude_age, !exclude_lt3_injections, !exclude_RVO, !exclude_DR_DMO, !exclude_no_va) %>% 
  select(-matches("exclude")) 


## Prepare time series of VA, OCT and fluid measurements for the selected eyes ##

# Just the VA history relating to the selected eyes
va_history <- va_raw %>% 
  inner_join(eye %>% 
               select(PatientID, EyeCode), by = c("PatientID", "EyeCode")) %>% 
  # Roll forward baseline VA measurements so that all series start at the index date 
  # (max roll forward is 2 weeks)
  mutate(months_since_index = if_else(baseline & months_since_index < 0, 0, months_since_index))


# Fluid history for each eye
fluid_history <- fluid_raw %>% 
  inner_join(eye %>% 
               select(PatientID, EyeCode), by = c("PatientID", "EyeCode"))


# OCT thickness history relating to the selected eyes
thickness_history <- oct_thickness_raw %>%  
  inner_join(eye %>% 
               select(PatientID, EyeCode), by = c("PatientID", "EyeCode")) %>% 
  # Calculate year of treatment for each eye (since index date)
  mutate(treatment_year = as.factor(floor(as.interval(index_date, ExamDate)/dyears()) +1 ))



## Calculate treatment intervals for the selected eyes

# Just the injections relating to selected eyes
injections <- injections_clean %>% 
  inner_join(eye %>% 
               select(PatientID, EyeCode, index_date), by = c("PatientID", "EyeCode")) %>% 
  # Calculate year of treatment for each eye (since index date)
  # Starts count at year 1
  mutate(treatment_year = as.factor(floor(as.interval(index_date, EncounterDate)/dyears()) +1))

# Note that some injections for selected eyes were in other clinic types e.g. cataract or glaucoma
# This may be because an incorrect clinic type was set at the beginning of the day or because of a genuine clinical need
# These injection dates are EXCLUDED from the NOA analysis but INCLUDED in the injection interval analysis
encounters_raw %>% 
  left_join(injections, by = "EncounterID") %>% 
  count(ClinicalCategoryDesc, EncounterTypeDesc, EncounterTypeCode, InjectedDrugDesc)


treatment_intervals <- injections %>%
  group_by(PatientID, EyeCode) %>% 
  mutate(next_encounter = lead(EncounterDate),
         treatment_interval_days = interval(EncounterDate, next_encounter)/ddays(), 
         treatment_interval_weeks = interval(EncounterDate, next_encounter)/dweeks(),
         # Indicator of short treatment interval 
         treatment_interval_short = treatment_interval_days <= 40) %>% 
  filter(!is.na(next_encounter)) %>% 
  # Timeline since the index date for each person
  # First and second year since treatment start concept not relevant because do not know the date of first injection (no clinic entry date recorded)
  mutate("Treatment interval" = if_else(treatment_interval_weeks > 12, ">12 weeks", "<=12 weeks"),
         # Calculate month and year of treatment for each eye (since index date)       
         treatment_months = (interval(index_date, EncounterDate)/dmonths()) + 1,
         treatment_year = as.factor(floor(as.interval(index_date, EncounterDate)/dyears()) + 1)) %>% 
  # Injection sequence (index injection = 1)
  mutate(injection_seq = row_number()) %>% 
  ungroup() %>% 
  mutate(sequence_id = abbreviate(paste(PatientID, EyeCode))) 

