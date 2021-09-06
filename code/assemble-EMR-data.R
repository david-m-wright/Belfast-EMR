# Assemble EMR dataset from Medisoft extract

library(tidyverse)
library(lubridate)
library(rprojroot)
library(eye)
library(TraMineR)
library(conflicted)
conflict_prefer("filter", "dplyr")

source(find_rstudio_root_file("code/EMR-helper-functions.R"))

## Load data ##

file_path <- find_rstudio_root_file("data/DMO data") 

# Load input files
emr_tables <- list.files(path = file_path, full.names = T) %>% 
  set_names(., str_remove(basename(.), "\\.txt")) %>%
  map(~read_delim(., delim = "|"))

# Description of the input files
emr_tables_desc <- read_csv(find_rstudio_root_file("study_design/Belfast-EMR-table-descriptions.csv")) %>% 
  mutate(across(`Table name`, str_replace, "^BEL", "DMO"))

# Pre-process records

# Patient details
patients <- emr_tables$DMOPatientDetails

# DMO diagnosis
diagnoses <- emr_tables$DMODiagnoses
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
dmo_encounters <- emr_tables$DMOEncounters

# Visual acuity
visual_acuity <- emr_tables$DMOVisualAcuity %>% 
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

visual_acuity %>% 
  count(PatientID, EyeCode, EncounterID) %>% 
  filter(n>1)

# Injections
injections_raw <- emr_tables$DMOInjections %>% 
  # Exclude non anti-VEGF injections
  mutate(exclude_not_antiVEGF = AntiVEGFInjection != 1) %>% 
  # Exclude a small number of duplicate encounters on the same date so only one injection per eye per date
  group_by(PatientID, EyeCode, EncounterDate, exclude_not_antiVEGF) %>% 
  summarise(row_number = row_number(), 
            EncounterID = EncounterID[row_number==1], .groups = "drop") %>% 
  mutate(exclude_duplicate_encounter = row_number != 1) %>% 
  select(-row_number) 

injections_clean <- injections_raw %>% 
  filter(!exclude_not_antiVEGF, !exclude_duplicate_encounter) %>% 
  
  # Add visual acuity at the selected encounter
  left_join(visual_acuity %>%
              select(PatientID, EyeCode, EncounterID, va_etdrs, va_logmar),
            by = c("PatientID", "EyeCode", "EncounterID")) 
      
## Assemble cohort (eye level) ##


# Injection history summary by eye
injection_summary_eye  <- injections_clean %>% 
  
  # Find date of first injection (IndexDate) and final injection
  # and number of injections for each eye
  group_by(PatientID, EyeCode) %>% 
  summarise(index_date = min(EncounterDate),
            final_injection_date = max(EncounterDate), 
            total_injections = n(),
            .groups = "drop") %>% 
  mutate(total_intervals = total_injections - 1) %>% 

  mutate(Injections = cut(total_injections, breaks = c(0,1,2,3, 12, max(total_injections)))) %>% 
  mutate(across(Injections, ~fct_relabel(., .fun = IntervalToInequality, unit = "n")))



# DMO diagnosis history of each eye
# Find most recent diagnosis preceding each injection
DMO_history_eye <- injections_clean %>% 
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
  left_join(injections_clean %>% 
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
  # filter(!exclude_age, !exclude_no_intervals)
  filter(!exclude_age)


# Just the injections relating to selected eyes
injections <- injections_clean %>% 
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
         # Month in which the interval started
          follow_up_months = interval(index_date, EncounterDate)/dmonths(),
         # Year of treatment in which the interval started 
         year_of_treatment = factor(floor(interval(index_date, EncounterDate)/dyears())+1),
         year_of_treatment = fct_collapse(year_of_treatment, ">=6" = c("6","7","8","9","10", "11"))) %>% 
  # Injection sequence (index injection = 1)
  mutate(injection_seq = row_number()) %>% 
  ungroup() %>% 
  mutate(sequence_id = abbreviate(paste(PatientID, EyeCode))) 

treatment_intervals %>% 
  count(follow_up_months, year_of_treatment) %>% print(n=100)

# Summarise treatment intervals by year of treatment
treatment_intervals_summary <-  treatment_intervals %>% 
  group_by(PatientID, EyeCode, year_of_treatment) %>% 
  summarise(any_gt12_weeks = any(`Treatment interval` == ">12 weeks"),
            n_gt12_weeks = sum(treatment_interval_weeks > 12),
            prop_gt12_weeks = n_gt12_weeks/length(treatment_interval_weeks), .groups = "drop")
  

# Find treatment naive patients

treatment_intervals_summary %>% 
  ggplot(aes(x = prop_gt12_weeks, fill =any_gt12_weeks)) +
  geom_histogram() +
  facet_wrap(~year_of_treatment)


# Proportion with any >12 weeks by year of treatment

# Baseline characteristics of those with any >12 weeks by year of treatment


# Would be useful to look at characteristics at start of year of treatment


clinical <- emr_tables$DMOClinicalFindings

co_pathology <- emr_tables$DMOCoPathology

co_pathology %>% 
  group_by(PatientID, EyeCode) %>% 
  count()

diabetic_diagnosis <- emr_tables$DMODiabeticDiagnosis

medical_history <- emr_tables$DMOMedicalHistory
  
# diabetic_diagnosis %>% select(-ExtractDate, -EncounterID) %>% count(DiabetesTypeDesc, AgeDiagnosed) %>% print(n=Inf)


# Assemble individual level summaries

indiv <- eye %>% 
  group_by(PatientID, Gender, index_age) %>%  
  summarise(Treated = if_else(n() == 2, "Bilateral", "Unilateral"), .groups = "drop")
  


# Years observed and years treated
eye %>% 
  summarise(across(c(years_observed, years_treated), sum))
  

# Check observation periods
# Note that for many patients the observation period is longer.
# This may be because they have now finished treatment.
# eye %>% 
#   ggplot(aes(x = years_observed, y = years_treated)) +
#   geom_point()
