# Assemble EMR dataset from Medisoft extract

library(tidyverse)
library(lubridate)
library(rprojroot)
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

# Injections
injections <- read_delim(paste0(file_path, "/DMOInjections.txt"), delim = "|") %>% 
  filter(AntiVEGFInjection == 1) %>% 
  # Remove a small number of duplicates
  group_by(PatientID, EyeCode, EncounterDate) %>% 
  summarise(row_number = row_number(), .groups = "drop") %>% 
  filter(row_number == 1) %>% 
  select(-row_number)
  

## Assemble cohort (eye level) ##

# Injection history summary by eye
# Find date of first injection (IndexDate) and final injection
# and number of injections for each eye
injection_summary_eye  <- injections %>% 
  group_by(PatientID, EyeCode) %>% 
  summarise(index_date = min(EncounterDate),
            final_injection_date = max(EncounterDate), 
            total_injections = n(),
            .groups = "drop") %>% 
  mutate(total_intervals = total_injections - 1)


# DMO diagnosis history of each eye
# Find most recent diagnosis preceding each injection
DMO_history_eye <- injections %>% 
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
  
  # Exclusion criteria
  mutate(
       # 18 or older on date of first injection (index date)
       exclude_age = index_age < 18 |is.na(index_age),
       
       # No DMO diagnosis within previous 180 days
       exclude_no_diagnosis = is.na(days_since_diag) | days_since_diag > 180,
       
       # No injections  after index date (so no treatment intervals can be calculated)
       exclude_no_intervals = total_intervals == 0)


# Apply the exclusion criteria
eye <- eye_raw %>% 
  filter(!exclude_age, !exclude_no_intervals)



# Assemble individual level summaries

indiv <- eye %>% 
  group_by(PatientID, Gender, index_age) %>%  
  summarise(Treated = if_else(n() == 2, "Bilateral", "Unilateral"), .groups = "drop")
  


