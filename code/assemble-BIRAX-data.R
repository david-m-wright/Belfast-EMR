# Assemble EMR dataset from Medisoft extract

library(tidyverse)
library(lubridate)
library(rprojroot)
library(eye)
# library(TraMineR)
library(conflicted)
conflict_prefer("filter", "dplyr")

source(find_rstudio_root_file("code/EMR-helper-functions.R"))

## Load data ##

file_path <- find_rstudio_root_file("data/BIRAX data") 


# Load input files
emr_table_list <-
  list.files(
    path = file_path,
    full.names = F,
    recursive = T,
    pattern = "\\.txt"
  ) %>% basename() %>%
  str_remove(pattern = "\\.txt") %>% 
  unique() %>% 
  as.list() %>% 
  set_names(., .)

emr_tables <- emr_table_list %>% 
    map( ~ read_delim(paste0(file_path, "/", ., ".txt"), delim = "|", guess_max = 10000))

# Description of the input files
emr_tables_desc <- read_csv(find_rstudio_root_file("data-dictionary/Belfast-EMR-table-descriptions.csv"))


# Pre-process records

# Clinical findings
clinical <- emr_tables$BELClinicalFindings

# CoPathology
# Output table for labelling which conditions affect visual fields and should be excluded from history or follow up
# emr_tables$BELCoPathology %>%
#   filter(CoPathologyDesc != "Glaucoma") %>%
#   count(CoPathologyDesc, CoPathologyCode) %>%
#   arrange(desc(n)) %>%
#   write_csv(file = find_rstudio_root_file("data-dictionary", "BEL-co-pathology-lookup-to-classify.csv"))

# Two columns, exclude history, exclude followup

copath <- emr_tables$BELCoPathology #%>%
#   # Classify co-pathologies as seriously affecting vision
#   inner_join(read_csv(
#     find_rstudio_root_file("data-dictionary", "GRIP-co-pathology-lookup.csv")
#   ),
#   by = c("CoPathologyDesc", "CoPathologyCode"))
# 
# # For each eye flag date of earliest serious co-pathology (i.e. affecting vision)
# copath_start <- copath %>%
#   filter(exclude_history | exclude_followup) %>%
#   group_by(PatientID, EyeCode, exclude_history, exclude_followup) %>%
#   summarise(copath_start = min(EncounterDate), .groups = "drop") %>%
#   rename(exclude_copath_history = exclude_history,
#          exclude_copath_followup = exclude_followup)
# 


# Diabetic Diagnosis
diabetes <- emr_tables$BELDiabeticDiagnosis %>%
  # Classify those with diabetes
  inner_join(read_csv(
    find_rstudio_root_file("data-dictionary", "Belfast-diabetes-lookup.csv")
  ),
  by = c("DiabetesTypeDesc", "DiabetesType"))

# Diabetes status
diab_status <- diabetes %>%
  filter(Diabetes) %>%
  group_by(PatientID, Diabetes) %>%
  # If there are multiple diabetes types recorded for one patient, take the most commonly recorded type
  summarise(
    DiabetesTypeDesc = MaxVotes(DiabetesTypeDesc),
    # Find youngest age at diagnosis
    age_diabetes_diagnosed = na_if(min(AgeDiagnosed, na.rm = T), Inf),
    .groups = "drop"
  )



# Diagnosis
# Excluding glaucoma diagnoses as these will be determined using visual fields
# emr_tables$GRIPDiagnoses %>%
#   filter(DiagnosisCategoryCode != "GLA") %>%
#   count(DiagnosisCategoryDesc, DiagnosisDescription) %>%
#   arrange(DiagnosisCategoryDesc, desc(n)) %>%
#   write_csv(file = find_rstudio_root_file("Data dictionary", "GRIP-co-pathology-lookup-to-classify.csv"))

diagnoses <- emr_tables$BELDiagnoses #%>%
#   inner_join(
#     read_csv(
#       find_rstudio_root_file("data-dictionary", "Belfast-diagnosis-lookup.csv")
#     ) %>%
#       select(-n),
#     by = c("DiagnosisCategoryDesc", "DiagnosisDescription")
#   )
# 
# # For each eye flag date of earliest serious co-pathology (i.e. affecting vision) excluding glaucoma
# # Glaucoma will be determined using visual fields
# diagnosis_start <- diagnoses %>%
#   filter(exclude_history | exclude_followup) %>%
#   group_by(PatientID, EyeCode, exclude_history, exclude_followup) %>%
#   summarise(diagnosis_start = min(DateofDiagnosis),
#             .groups = "drop") %>%
#   rename(exclude_diagnosis_history = exclude_history,
#          exclude_diagnosis_followup = exclude_followup)


# DR assessment
dr_assessment <- emr_tables$BELDRAssessment

# DR grading
dr_grading <- emr_tables$BELDRGrading

# Encounters
# These will be important for health economics evaluation
# Relevant encounters for patient history are visual fields, dates of diagnosis of any of the co-pathologies affecting visual fields.
encounters <- emr_tables$BELEncounters

# Investigations
investigations <- emr_tables$BELInvestigations

# IOP
iop <- emr_tables$BELIOP %>%
  rename(IOP = Value)


# Medical history
med_history  <- emr_tables$BELMedicalHistory

# Medication
medication <- emr_tables$BELMedication

# Post operative complications
post_op <- emr_tables$BELPostOperativeComplications

#  Load patient details
patients <- emr_tables$BELPatientDetails %>%
  mutate(across(c(Gender, EthnicDesc),  ~fct_explicit_na(as.factor(.), na_level = "Missing")))

# Surgery
surgery <- emr_tables$BELSurgery

# Surgery indications
surgery_ind <- emr_tables$BELSurgeryIndications

# Visual acuity
visual_acuity <- emr_tables$BELVisualAcuity %>% 
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
injections_raw <- emr_tables$BELInjections %>% 
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
  select(-exclude_not_antiVEGF, -exclude_duplicate_encounter) %>%   
  # Add visual acuity at the selected encounter
  left_join(visual_acuity %>%
              select(PatientID, EyeCode, EncounterID, va_etdrs, va_logmar),
            by = c("PatientID", "EyeCode", "EncounterID")) 

# OCT thickness maps
oct_thickness <- emr_tables$OCT_ThicknessMetrics



## Assemble cohort (eye level) ##

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
  
  # # Check DMO diagnosis within previous 180 days
  # left_join(DMO_history_eye,
  #           by = c("PatientID" = "PatientID", "EyeCode" = "EyeCode", "index_date" = "EncounterDate")) %>% 
  
  # Add visual acuity at baseline
  left_join(injections_clean %>% 
              rename("ETDRS letters" = va_etdrs),
            by = c("PatientID","EyeCode", "index_date" = "EncounterDate")) %>% 
  
  # Exclusion criteria
  mutate(
    # 18 or older on date of first injection (index date)
    exclude_age = index_age < 18 |is.na(index_age),
    
    # # No DMO diagnosis within previous 180 days
    # exclude_no_diagnosis = is.na(days_since_diag) | days_since_diag > 180,
    
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
         follow_up_months = interval(index_date, EncounterDate)/dmonths()) %>% 
  # Injection sequence (index injection = 1)
  mutate(injection_seq = row_number()) %>% 
  ungroup() %>% 
  mutate(sequence_id = abbreviate(paste(PatientID, EyeCode))) 


# Find patients with an AMD diagnosis
# AMD_patients <- 


# Find treatment naive patients
  

# 
# # What about VA measurements taken on the same day?
# # Produces sawtooth pattern that should be removed
# # Which measurements to prioritise?
# visual_acuity %>% 
#   slice(1000) %>% 
#   select(PatientID) %>% 
#   inner_join(visual_acuity, by = "PatientID") %>%
#   filter(RecordedNotation == "LETTERSCORE") #%>% 
# ggplot(aes(x = EncounterDateTime, y = va_logmar, col = EyeCode)) +
#   geom_line() +
#   scale_y_reverse() +
#   theme_light()
# 
# 
# 
# visual_acuity %>% 
#   slice(3) %>% 
#   select(PatientID) %>% 
#   inner_join(visual_acuity, by = "PatientID") %>% 
#   filter(EyeCode=="L") %>% 
#   # select(EncounterDate, va_logmar) %>%
#   select(EncounterDateTime, RecordedNotation, RecordedNotationBestCorrected, snellen_to_logmar, etdrs_to_logmar, va_logmar) %>% 
#   arrange(EncounterDateTime) %>% 
#   print(n=Inf)
# 
# # Sequences of VA are long
# visual_acuity %>% 
#   count(PatientID) %>%
#   pull(n) %>% 
#   table() %>% 
#   plot()
#   table()
#   ggplot(aes(x = n)) +
#   geom_histogram()
#   
  
  