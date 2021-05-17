# Assemble EMR dataset from Medisoft extract

library(tidyverse)
library(lubridate)
library(rprojroot)
library(eye)
library(TraMineR)
# Conflicted?


source(find_rstudio_root_file("code/EMR-helper-functions.R"))

## Load data ##

file_path <- find_rstudio_root_file("data/DMO data") 

#  Load patient details
patients <- read_delim(paste0(file_path, "/DMOPatientDetails.txt"), delim = "|")

diagnoses <- read_delim(paste0(file_path, "/DMODiagnoses.txt"), delim = "|")

# DMO diagnosis
diagnoses %>% 
  count(DiagnosisDescription, str_detect(DiagnosisICD10Code, "H35.81")) %>% print(n=Inf)

diagnoses %>% 
  mutate(oed = DiagnosisICD10Code == "H36.0 + E14.3"
         & str_detect(DiagnosisDescription, "oedema|maculopathy")) %>% 
  # ) %>% 
  count(DiagnosisICD10Code, DiagnosisDescription, oed) %>% 
  arrange(n) %>% 
  print(n=Inf)

# DMO encounters
dmo_encounters <- read_delim(paste0(file_path, "/DMOEncounters.txt"), delim = "|")

# Visual acuity
visual_acuity <- read_delim(paste0(file_path, "/DMOVisualAcuity.txt"), delim = "|") %>% 
  # Converting visual acuity to ETDRS 
  mutate(cleaned_etdrs = if_else(str_detect(RecordedNotation, "LETTERSCORE"), 
                                 as.character(va(RecordedNotationBestMeasure, from = "etdrs", to = "etdrs")), as.character(NA)),
         snellen_to_etdrs = case_when(str_detect(RecordedNotation, "SNELLEN(FEET|METRES)") ~ 
                                        as.character(va(RecordedNotationBestMeasure, from = "snellen", to = "etdrs")),
                                      str_detect(RecordedNotation, "SNELLENFRACTION") ~ 
                                        as.character(va(RecordedNotationBestMeasure, from = "snellendec", to = "etdrs"))),
         # Note that the 2DP logmar is rounded to 1DP before converion to ETDRS
         logmar_to_etdrs = as.character(va(coalesce(Logmar2DBestMeasure, Logmar1DBestMeasure), from = "logmar", to = "etdrs")),
         va_etdrs = as.numeric(coalesce(cleaned_etdrs, snellen_to_etdrs, logmar_to_etdrs))) %>% 
  
  mutate(cleaned_logmar = as.character(va(coalesce(Logmar2DBestMeasure, Logmar1DBestMeasure), from = "logmar", to = "logmar")),
         etdrs_to_logmar = if_else(str_detect(RecordedNotation, "LETTERSCORE"), 
                                   as.character(va(RecordedNotationBestMeasure, from = "etdrs", to = "logmar")), as.character(NA)),
         snellen_to_logmar = case_when(str_detect(RecordedNotation, "SNELLEN(FEET|METRES)") ~ 
                                         as.character(va(RecordedNotationBestMeasure, from = "snellen", to = "logmar")),
                                       str_detect(RecordedNotation, "SNELLENFRACTION") ~ 
                                         as.character(va(RecordedNotationBestMeasure, from = "snellendec", to = "logmar"))),
         va_logmar = as.numeric(coalesce(cleaned_logmar, snellen_to_logmar, etdrs_to_logmar))) 

# Injections
injections_raw <- read_delim(paste0(file_path, "/DMOInjections.txt"), delim = "|") %>% 
  filter(AntiVEGFInjection == 1) %>% 
  # Remove a small number of duplicate encounters on the same date
  group_by(PatientID, EyeCode, EncounterDate) %>% 
  summarise(row_number = row_number(), 
            EncounterID = EncounterID[row_number==1], .groups = "drop") %>% 
  filter(row_number == 1) %>% 
  select(-row_number) %>% 

  # Added visual acuity at the same encounter
  left_join(visual_acuity %>% 
            select(PatientID, EyeCode, EncounterID, va_etdrs, va_logmar),
            by = c("PatientID", "EyeCode", "EncounterID")) %>% 
  
  # Remove some more duplicate encounters so only one injection per eye per date
  group_by(PatientID, EyeCode, EncounterDate, va_etdrs, va_logmar) %>% 
  summarise(row_number = row_number(), .groups = "drop") %>% 
  filter(row_number == 1) %>% 
  select(-row_number) 
  



      
## Assemble cohort (eye level) ##

# Injection history summary by eye
# Find date of first injection (IndexDate) and final injection
# and number of injections for each eye
injection_summary_eye  <- injections_raw %>% 
  group_by(PatientID, EyeCode) %>% 
  summarise(index_date = min(EncounterDate),
            final_injection_date = max(EncounterDate), 
            total_injections = n(),
            .groups = "drop") %>% 
  mutate(total_intervals = total_injections - 1)


# DMO diagnosis history of each eye
# Find most recent diagnosis preceding each injection
DMO_history_eye <- injections_raw %>% 
  left_join(diagnoses %>% 
               filter(DiagnosisICD10Code == "H36.0 + E14.3") %>% 
               distinct(PatientID, DateofDiagnosis, EyeCode),
                        by = c("PatientID", "EyeCode")) %>% 
  mutate(DateofDiagnosis = if_else(DateofDiagnosis > EncounterDate, as.POSIXct(NA), DateofDiagnosis)) %>% 
  group_by(PatientID, EyeCode, EncounterDate) %>% 
  summarise(DateofDiagnosis = max(DateofDiagnosis), .groups = "drop") %>% 
  mutate(days_since_diag = interval(DateofDiagnosis, EncounterDate)/ddays(1))
  
  

# Assemble eye level dataset (single row per eye) 

eye_raw <- patients %>% 
  
  # Calculate age at first injection (index date)
  inner_join(injection_summary_eye,
            by = "PatientID") %>% 
  mutate(index_age =  as.interval(PerturbedDateofBirth, index_date)/dyears()) %>% 
  # Calculate years observed (index to current age/death)
  mutate(years_observed = coalesce(PerturbedCurrentAge, PerturbedAgeatDeath) - index_age) %>% 
  # Calculate years treated (index to final injection)
  mutate(years_treated = interval(index_date, final_injection_date)/dyears()) %>% 
  
  # Check DMO diagnosis within previous 180 days
  left_join(DMO_history_eye,
            by = c("PatientID" = "PatientID", "EyeCode" = "EyeCode", "index_date" = "EncounterDate")) %>% 
  
  # Add visual acuity at baseline
  left_join(injections_raw %>% 
              rename("ETDRS letters" = va_etdrs),
            by = c("PatientID","EyeCode", "index_date" = "EncounterDate")) %>% 
  
  # Exclusion criteria
  mutate(
       # 18 or older on date of first injection (index date)
       exclude_age = index_age < 18 |is.na(index_age),
       
       # No DMO diagnosis within previous 180 days
       exclude_no_diagnosis = is.na(days_since_diag) | days_since_diag > 180,
       
       # No injections  after index date (so no treatment intervals can be calculated)
       exclude_no_intervals = total_intervals == 0) %>% 
  
  # Derived variables
  mutate("ETDRS category" = fct_explicit_na(cut(`ETDRS letters`, 
                                breaks = c(0, 32, 73, 100), 
                                labels = c("<33 letters", "33-73 letters", ">73 letters"), include.lowest = T), 
                                na_level = "Missing"))

# Apply the exclusion criteria
eye <- eye_raw %>% 
  filter(!exclude_age, !exclude_no_intervals)

# Just the injections relating to selected eyes
injections <- injections_raw %>% 
  inner_join(eye %>% 
               select(PatientID, EyeCode, index_date), by = c("PatientID", "EyeCode")) 

## Calculate treatment intervals for the cohort ##

treatment_intervals <- injections %>%
  group_by(PatientID, EyeCode) %>% 
  mutate(next_encounter = lead(EncounterDate),
         treatment_interval_weeks = interval(EncounterDate, next_encounter)/dweeks()) %>% 
  filter(!is.na(next_encounter)) %>% 
  # Timeline since the index date for each person
  # First and second year since treatment start concept not relevant because do not know the date of first injection (no clinic entry date recorded)
    mutate("Treatment interval" = if_else(treatment_interval_weeks > 12, ">12 weeks", "<=12 weeks"),
         follow_up_months = interval(index_date, EncounterDate)/dmonths()) %>% 
  # Injection sequence (index injection = 1)
  mutate(injection_seq = row_number()) %>% 
  ungroup() %>% 
  mutate(sequence_id = abbreviate(paste(PatientID, EyeCode))) 


# Find treatment naive patients
  

clinical <- read_delim(paste0(file_path, "/DMOClinicalFindings.txt"), delim = "|")

co_pathology <- read_delim(paste0(file_path, "/DMOCoPathology.txt"), delim = "|")

diabetic_diagnosis <- read_delim(paste0(file_path, "/DMODiabeticDiagnosis.txt"), delim = "|")

medical_history <- read_delim(paste0(file_path, "/DMOMedicalHistory.txt"), delim = "|")
  
# diabetic_diagnosis %>% select(-ExtractDate, -EncounterID) %>% count(DiabetesTypeDesc, AgeDiagnosed) %>% print(n=Inf)


# Assemble individual level summaries

indiv <- eye %>% 
  group_by(PatientID, Gender, index_age) %>%  
  summarise(Treated = if_else(n() == 2, "Bilateral", "Unilateral"), .groups = "drop")
  




