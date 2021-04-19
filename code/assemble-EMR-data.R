# Assemble EMR dataset from Medisoft extract

library(tidyverse)
library(lubridate)
library(rprojroot)
# Conflicted?


source(find_rstudio_root_file("code/EMR-helper-functions.R"))

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
injections <- read_delim(paste0(file_path, "/DMOInjections.txt"), delim = "|") 

# Find date of first injection (IndexDate)
index_date <- injections %>% 
  filter(AntiVEGFInjection == 1) %>% 
  group_by(PatientID) %>% 
  summarise(IndexDate = min(EncounterDate))
  

# DMO diagnosis history of each patient
DMO_history <- index_date %>% 
  
  # Find dates of diagnosis with DMO (multiple records per eye)
  left_join(diagnoses %>% 
              filter(DiagnosisICD10Code == "H36.0 + E14.3") %>% 
              distinct(PatientID, DateofDiagnosis, EyeCode),
            by = "PatientID") %>% 
  
  # Was DMO diagnosis within previous 180 days of index date?
  mutate(Previous180 = interval(IndexDate - ddays(x=180), IndexDate),
         DMODiagnosis180 = DateofDiagnosis %within% Previous180)

  


# Injection history by eye
injection_history <- injections %>% 
  inner_join(index_date, by = "PatientID") %>% 
  filter(EncounterDate >= IndexDate) %>% 
  select(PatientID, EncounterDate, EyeCode) %>% 
  pivot_wider(names_from = EyeCode, values_from = EyeCode, id_cols = c(PatientID, EncounterDate), values_fn = length, values_fill = 0)




# eye_treated <- injections %>% 
#   count(PatientID, EncounterDate) %>% 
#   filter(n <=2) %>% 
#   mutate(Treated = if_else(n == 1, "Unilateral", "Bilateral"))



# Exclusion due to other conditions in previous 180 days ****

# Assemble individual level dataset            
indiv_raw <- patients %>% 
  
  # Calculate age at first injection (index date)
  left_join(index_date,
            by = "PatientID") %>% 
  mutate(IndexAge =  as.interval(PerturbedDateofBirth, IndexDate)/dyears()) %>% 
  # Calculate years observed (index to current age/death)
  mutate(YearsObserved = coalesce(PerturbedCurrentAge, PerturbedAgeatDeath) - IndexAge) %>% 
  
  # Check DMO diagnosis within previous 180 days
  left_join(DMO_history %>% 
      group_by(PatientID) %>% 
      summarise(DMODiagnosis180 = any(DMODiagnosis180)) ,
            by = "PatientID") %>% 
  
  # Laterality treated at index date
  left_join(injection_history %>% 
              mutate(Treated = if_else(L+R ==1, "Unilateral", "Bilateral")) %>% 
              select(PatientID, EncounterDate, Treated),
            by = c("PatientID" = "PatientID", "IndexDate" = "EncounterDate")) %>% 
  
  # Number of injections in each eye (including at index date)
  left_join(injection_history %>% 
              group_by(PatientID) %>% 
              summarise(across(c(L, R), sum)),
            by = "PatientID"
  ) %>% 
  
  
  # Calculate age in years at first encounter for DMO
  # left_join(dmo_encounters %>% 
  #             group_by(PatientID) %>% 
  #             summarise(FirstEncounter = min(EncounterDate)),
  #           by = "PatientID") %>% 
  # mutate(ToFirstEncounter = as.interval(PerturbedDateofBirth, FirstEncounter),
  #         BaselineAge = ToFirstEncounter/dyears()) %>% 
  
  
  # Exclusion criteria
  # No injections recorded
  mutate(exclude_injections = is.na(IndexDate),
      
         # 18 or older on date of first injection (index date)
         exclude_age = IndexAge < 18 |is.na(IndexAge),
         
         # No DMO diagnosis within previous 180 days
         exclude_no_diagnosis = !DMODiagnosis180 | is.na(DMODiagnosis180),
         
         # No injections in either eye after index date (so no intervals can be calculated)
         exclude_no_intervals = L <= 1 & R <= 1)
         
         # First encounter after the extraction date (error induced by perturbation?)
         # exclude_future = FirstEncounter > ExtractDate,
         # exclude_years_observed = YearsObserved < 0)


indiv_raw %>% count(exclude_no_intervals, Treated)
indiv_raw %>% 
  count(exclude_injections, exclude_age )

summary(indiv_raw$IndexAge)


# Individual level dataset 
indiv <- indiv_raw %>% 
  
  # Apply exclusion criteria
  filter(!exclude_injections,
         !exclude_age,
         !exclude_no_diagnosis,
         !exclude_no_intervals)

# Number of eyes
n_eyes <- nrow(filter(indiv, R > 1)) + nrow(filter(indiv, L > 1))
    
# Number of intervals 
n_intervals <- indiv %>% 
  mutate(across(c(L, R), ~abs(.-1))) %>% 
  summarise(across(c(L, R), sum)) %>% 
  mutate(total = L + R) %>% 
  pull(total)

# indiv %>% 
# ggplot(aes(x = YearsObserved)) +
#     geom_histogram(bins=30)
#   mutate(BaselineAge = FirstEncounter - PerturbedDateofBirth)


    # summarise(across(where(is.numeric), SummariseSD)) %>% 
  