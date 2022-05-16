# Assemble EMR dataset from Medisoft extract

library(tidyverse)
library(lubridate)
# data.table for speed in lead/lag operations
library(data.table)
library(rprojroot)
library(eye)
# library(TraMineR)
library(conflicted)
conflict_prefer("filter", "dplyr")

source(find_rstudio_root_file("code/EMR-helper-functions.R"))

## Load data ##

# Derived variable labels
var_desc <- read_csv(find_rstudio_root_file("data-dictionary/Belfast-EMR-derived-variables.csv"))


file_path <- find_rstudio_root_file("data/BIRAX data") 

# Characteristics of input files
emr_table_list <- 
list.files(
  path = file_path,
  full.names = F,
  recursive = T,
  pattern = "\\.txt"
) %>% basename() %>%
  str_remove(pattern = "\\.txt") %>% 
  unique() %>% 
  enframe(name = "Table number", value = "Table name") %>% 
  mutate(Rows = map_dbl(`Table name`, ~nrow(fread(paste0(file_path, "/", ., ".txt"), select = 1L))),
         Columns = map_dbl(`Table name`, ~ncol(fread(paste0(file_path, "/", ., ".txt"), sep = "|", nrows = 1L)))) %>% 
  inner_join(read_csv(find_rstudio_root_file("data-dictionary/Belfast-EMR-table-descriptions.csv")), by = "Table name")

# Description of the input files
# emr_tables_desc <- read_csv(find_rstudio_root_file("data-dictionary/Belfast-EMR-table-descriptions.csv"))

# Pre-process records

# Clinical findings
clinical <- fread(file.path(file_path, "BELClinicalFindings.txt"))

# CoPathology
# Output table for labelling which conditions affect visual fields and should be excluded from history or follow up
# emr_tables$BELCoPathology %>%
#   filter(CoPathologyDesc != "Glaucoma") %>%
#   count(CoPathologyDesc, CoPathologyCode) %>%
#   arrange(desc(n)) %>%
#   write_csv(file = find_rstudio_root_file("data-dictionary", "BEL-co-pathology-lookup-to-classify.csv"))

# Two columns, exclude history, exclude followup

copath <- fread(file.path(file_path, "BELCoPathology.txt"))
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
diabetes <- fread(file.path(file_path, "BELDiabeticDiagnosis.txt")) %>%
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

diagnoses <- fread(file.path(file_path, "BELDiagnoses.txt")) #%>%
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
dr_assessment <- fread(file.path(file_path, "BELDRAssessment.txt")) #%>%

# DR grading
dr_grading <- fread(file.path(file_path, "BELDRGrading.txt")) 

# Encounters
# These will be important for health economics evaluation
# Relevant encounters for patient history are visual fields, dates of diagnosis of any of the co-pathologies affecting visual fields.
encounters_raw <- fread(file.path(file_path, "BELEncounters.txt"))


# Exclude encounters that originated from other clinics (i.e. glaucoma service)
encounters <- encounters_raw %>% 
  filter(ClinicalCategoryDesc == "Medical retina",
         # Retain only the three major types of appointment in medical retina (nurse, doctor, operation/injection)
         EncounterTypeCode %in% c("SCLV", "NCLV", "OPER"))


# Investigations
investigations <- fread(file.path(file_path, "BELInvestigations.txt")) 

# IOP
iop <- fread(file.path(file_path, "BELIOP.txt")) %>% 
  rename(IOP = Value)


# Medical history
med_history  <- fread(file.path(file_path, "BELMedicalHistory.txt")) 

# Medication
medication <- fread(file.path(file_path, "BELMedication.txt")) 

# Post operative complications
post_op <- fread(file.path(file_path, "BELPostOperativeComplications.txt"))

#  Load patient details
patients <- fread(file.path(file_path, "BELPatientDetails.txt")) %>% 
  mutate(across(c(Gender, EthnicDesc),  ~fct_explicit_na(as.factor(.), na_level = "Missing")))

# Surgery
surgery <- fread(file.path(file_path, "BELSurgery.txt"))

# Surgery indications
surgery_ind <- fread(file.path(file_path, "BELSurgeryIndications.txt")) 

# Visual acuity
visual_acuity <- fread(file.path(file_path, "BELVisualAcuity.txt")) %>% 
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
         va_logmar = as.numeric(coalesce(cleaned_logmar, snellen_to_logmar, etdrs_to_logmar))) %>% 

  # Categorise VA into functionally relevant categories
  mutate(va_category_snellen = fct_explicit_na(cut(va_logmar, 
                  breaks = c(min(va_logmar, na.rm=TRUE), 0.29, 0.59, 0.99, max(va_logmar, na.rm=TRUE)), 
                  labels = c("Good (VA>6/12)", "Moderate (6/24 < VA <= 6/12)", "Partially sighted (6/60 < VA <= 6/24)", "Blind (VA <= 6/60)"), 
                  include.lowest = TRUE),
                  na_level = "Missing"), 
         va_category_etdrs = fct_explicit_na(cut(va_etdrs,
                                              breaks = c(0, 32, 73, 100),
                                              labels = c("<33 letters", "33-73 letters", ">73 letters"), include.lowest = T),
                                          na_level = "Missing")) 

# Injections
injections_raw <- fread(file.path(file_path, "BELInjections.txt")) %>% 
  # Exclude non anti-VEGF injections
  mutate(exclude_not_antiVEGF = AntiVEGFInjection != 1) %>% 
  # Exclude a small number of duplicate encounters on the same date so only one injection per eye per date
  group_by(PatientID, EyeCode, EncounterDate, InjectedDrugDesc, exclude_not_antiVEGF) %>% 
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
oct_thickness <- fread(file.path(file_path, "OCT_ThicknessMetrics.txt")) 



## Assemble cohort (eye level) ##


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
  slice_min(order_by = "DateOfDiagnosis_AMD", with_ties = FALSE) %>% 
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
  slice_min(order_by = "DateOfDiagnosis_RVO", with_ties = FALSE) %>% 
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
  slice_min(order_by = "DateOfDiagnosis_DR", with_ties = FALSE) %>% 
  ungroup()



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




### Take snapshots of VA at various stages of follow-up ###

# For each visual acuity measurement, calculate the time since the first injection (index_date)
# 72s
# tic()
va_long <- injection_summary_eye %>% 
  as.data.table() %>% 
  # Join injection_summary_eye and visual_acuity tables
  .[visual_acuity, .(PatientID, EyeCode, index_date, EncounterDate, va_logmar, va_etdrs, va_category_snellen), on = .(PatientID, EyeCode)] %>% 
  # Drop missing VA measurements
  .[!is.na(va_logmar)] %>% 
  # Retain VA measurements from up to one month prior to the index date (to allow for small mismatches in initial VA measurement and injection date)
  .[EncounterDate >= (index_date - dmonths(1))] %>% 
  .[, months_since_index:=as.interval(index_date, EncounterDate)/dmonths()] %>% 
  # Find the best VA recorded by date for each eye
  .[, .SD[which.min(va_logmar)], by = c("PatientID", "EyeCode", "months_since_index")]
# toc()




## Find closest VA measurement to each of the specified snapshot times, within set tolerances

# Follow-up times in months at which the VA status is to be taken, along with acceptable tolerances
snapshots <- data.table(months_since_index = c(0, 12, 24, 36, 48, 60, 84, 120), 
                        follow_up = c(0, 12, 24, 36, 48, 60, 84, 120), 
                        tolerance = c(1, 2, 2, 2, 2, 2,2,2),
                        key = "months_since_index") %>% 
  .[, `:=` (lower_lim = follow_up - tolerance,
            upper_lim = follow_up + tolerance,
            follow_up_year = as.factor(follow_up/12))]


# Set keys to join on (so the .on() specification is not needed in the join)
setkey(va_long, months_since_index)

# Find the nearest candidate VA measurement for each follow-up snapshot (snapshots are on injection time timescale)
# The first argument inside the square brackets gives the number of rows of the output (the index i)
# i.e. the lookup table is outside the brackets
va_history_raw <- snapshots[va_long, roll = "nearest" ] %>% 
  unique() %>% 
  .[, within_tolerance := months_since_index >= lower_lim & months_since_index <= upper_lim] %>% 
  # Find the nearest injection time to the desired follow-up snapshot (reverse of above)
  .[, nearest_to_follow_up := min(abs(months_since_index-follow_up)), by = c("PatientID", "EyeCode", "follow_up")] %>% 
  # Retain all measurements but mark those that can be used for a snapshot
  # Break  ties of equal distance from the snapshot by choosing the earlier VA measurement (rolling back)
  .[, nearest_within_tolerance := abs(months_since_index-follow_up) == nearest_to_follow_up & within_tolerance] %>% 
  .[, snapshot := abs(months_since_index-follow_up) == nearest_to_follow_up & within_tolerance & .I == min(.I), by = c("PatientID", "EyeCode", "follow_up", "nearest_within_tolerance")]

# Retain only those with a VA measurement close enough to the index date to be called a baseline
# Inner join to self
va_history_clean <- va_history_raw[va_history_raw[snapshot == T & follow_up == 0, .(PatientID, EyeCode)], on = .(PatientID, EyeCode)] %>% 
  # Calculate VA change from baseline snapshot (VA measurement closest to index date)
  setkeyv(., c("PatientID", "EyeCode", "EncounterDate")) %>% 
  .[,  `:=` (va_change_logmar = va_logmar - va_logmar[(snapshot & follow_up == 0)],
             va_change_etdrs = va_etdrs - va_etdrs[(snapshot & follow_up == 0)]),
    , by=.(PatientID, EyeCode)] %>% 
# .[, `:=` (va_change_gteq10 = factor(if_else(va_change_logmar <= -0.2, "Improvement", "No improvement"), levels = c("No improvement", "Improvement")),
#           va_change_gteq15 = factor(if_else(va_change_logmar <= -0.3, "Improvement", "No improvement"), levels = c("No improvement", "Improvement")))] %>% 
# VA change into three categories
mutate(va_change_lines = cut(va_change_logmar, 
                             breaks = c(min(va_change_logmar), -0.2, 0.19, max(va_change_logmar)), 
                             labels = c("gained >=2 lines", "< 2 lines change", "lost >= 2 lines"), include.lowest = TRUE))


# Assemble eye level dataset (single row per eye) 

eye_raw <- patients %>% 
  
  # Add diagnoses for AMD, RVO and DMO/DR
  left_join(AMD_patients, by = "PatientID") %>% 
  left_join(RVO_patients, by = c("PatientID", "EyeCode")) %>% 
  left_join(DMO_DR_patients, by = c("PatientID", "EyeCode")) %>% 
  # Those with a valid VA history
  left_join(va_history_clean %>% 
              distinct(PatientID, EyeCode, exclude_no_va = FALSE), by = c("PatientID", "EyeCode")) %>% 
  mutate(across(exclude_no_va, ~if_else(is.na(.), TRUE, .))) %>% 
  
  # Calculate age at first injection (index date)
  left_join(injection_summary_eye, by = c("PatientID", "EyeCode")) %>% 
  mutate(index_age =  as.interval(PerturbedDateofBirth, index_date)/dyears()) %>% 
  # Calculate years observed (index to current age/death)
  mutate(years_observed = coalesce(PerturbedCurrentAge, PerturbedAgeatDeath) - index_age) %>% 
  # Calculate years treated (index to final injection)
  mutate(years_treated = interval(index_date, final_injection_date)/dyears()) %>% 
  
  # Add visual acuity at baseline
  left_join(injections_clean %>% 
              rename("ETDRS letters" = va_etdrs) %>% 
              select(-EncounterID),
            by = c("PatientID","EyeCode", "index_date" = "EncounterDate")) %>% 
  
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
    exclude_no_intervals = total_intervals == 0) 


# Apply the exclusion criteria
eye <- eye_raw %>% 
  filter(!exclude_no_AMD, !exclude_age, !exclude_lt3_injections, !exclude_RVO, !exclude_DR_DMO, !exclude_no_va)
  

# Just the injections relating to selected eyes
injections <- injections_clean %>% 
  inner_join(eye %>% 
               select(PatientID, EyeCode, index_date), by = c("PatientID", "EyeCode")) %>% 
  # Calculate year of treatment for each eye (since index date)
  mutate(treatment_year = as.factor(floor(as.interval(index_date, EncounterDate)/dyears()) +1))

# Note that some injections for selected eyes were in other clinic types e.g. cataract or glaucoma
# This may be because an incorrect clinic type was set at the beginning of the day or because of a genuine clinical need
# These injection dates are EXCLUDED from the NOA analysis but INCLUDED in the injection interval analysis
encounters_raw %>% 
  left_join(injections, by = "EncounterID") %>% 
  count(ClinicalCategoryDesc, EncounterTypeDesc, EncounterTypeCode, InjectedDrugDesc)


# Just the VA history relating to the selected eyes
va_history <- va_history_clean %>% 
  inner_join(eye %>% 
               select(PatientID, EyeCode), by = c("PatientID", "EyeCode"))

# OCT visit history relating to the selected eyes
oct_history <- oct_thickness %>% 
  distinct(PatientID, EyeCode, ExamDate) %>% 
  inner_join(eye %>% 
               select(PatientID, EyeCode, index_date), by = c("PatientID", "EyeCode")) %>% 
  # Calculate year of treatment for each eye (since index date)
  mutate(treatment_year = as.factor(floor(as.interval(index_date, ExamDate)/dyears()) +1 ))

# Note many OCTs prior to treatment initiation.
# oct_history %>% count(treatment_year)


## Calculate treatment intervals for the selected eyes

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

# Metrics to indicate stability of patient based on treatment interval sequence
treatment_intervals %>% 
  select(sequence_id, EncounterDate, treatment_interval_weeks, treatment_months)

# Define patient as stable in a given year if 
# No treatment intervals <= 40 days?
# What about first year with loading doses - would classify all patients as unstable
treatment_intervals %>% 
  group_by(sequence_id) %>% 
  mutate(short_interval = treatment_interval_days <= 40,
    next_interval_short = dplyr::lead(treatment_interval_days <= 30)) %>% 
  group_by(sequence_id, treatment_year) %>% 
  mutate(longest_interval = max(treatment_interval_weeks),
         ratio_present_to_last = treatment_interval_weeks/dplyr::lag(treatment_interval_weeks)) %>% 
  select(sequence_id, treatment_interval_days, short_interval)
  #select(sequence_id, EncounterDate, treatment_interval_days, treatment_months, next_interval_short, longest_interval, ratio_present_to_last)



# NOA data dictionary
noa_dictionary <- fread(find_rstudio_root_file("data-dictionary/NOA_dictionary.csv")) 

# Progress logs from the NOA
latest_noa_log <- list.files(file.path(file_path, "NOA-output"), pattern = "[0-9]{2}-[0-9]{2}-[0-9]{4}") %>% 
  str_extract("[0-9]{2}-[0-9]{2}-[0-9]{4}") %>% 
  dmy() %>% 
  max()

# Latest log is extracted into the 'Logs' folder
noa_log <- fread(file.path(file_path, "NOA-output", "Logs", "Report.csv")) %>% 
  mutate(LogDate = latest_noa_log) %>% 
  relocate(LogDate)

# List of unsuccessful scans and reasons
noa_unsuccessful <- fread(file.path(file_path, "NOA-output", "Logs", "Unsuccessful_scans.csv")) %>% 
  separate(FileName, into = c("PatientID", "EyeCode", "OCTDate", "OCTSequence"), sep = "_", remove = FALSE) %>% 
  mutate(across(EyeCode, ~if_else(.== "OD", "R", "L")),
         across(OCTDate, ymd),
         across(OCTSequence, ~as.numeric(if_else(is.na(OCTSequence), "0", OCTSequence))))  


# Output from the Notal Ophthalmic Analyzer (NOA)
noa_raw <- fread(file.path(file_path, "NOA-output", "Results.csv"),
             na.strings = "nan") %>% 
  separate(FileName, into = c("PatientID", "EyeCode", "OCTDate", "OCTSequence"), sep = "_", remove = FALSE) %>% 
  mutate(across(EyeCode, ~if_else(.== "OD", "R", "L")),
         across(OCTDate, ymd),
         across(OCTSequence, ~as.numeric(if_else(is.na(OCTSequence), "0", OCTSequence))),
         across(matches("Bscan_|Section_"), ~na_if(., -1000)),
         across(c(IRF, 
                  SRF, 
                  Fluid, 
                  `Bscan_1st_highest_fluid[nl]`,
                  `Bscan_2nd_highest_fluid[nl]`,
                  `Bscan_3rd_highest_fluid[nl]`,
                  ERM,
                  RPE_irregularities,
                  PED), as.factor))


# Where there were multiple scans, select the best quality scan for each eye for each visit date
# Selected scans were eligible where possible
# Fewest missing values (indicator of neural network processing),
# Highest number of frames in the AVI (B-scans)
# If tied, the final scan taken on that date is selected.
noa <- noa_raw %>% 
  # filter(`Analysis eligibility` == 1) %>% 
  mutate(OCTMissing = rowSums(across(everything(), is.na))) %>% 
  group_by(PatientID, EyeCode, OCTDate) %>% 
  slice_max(`Analysis eligibility`, n = 1, with_ties = TRUE) %>% 
  slice_min(OCTMissing, n = 1, with_ties = TRUE) %>% 
  slice_max(FramesNum, n = 1, with_ties = TRUE) %>% 
  slice_max(OCTSequence, n = 1, with_ties = FALSE) %>% 
  ungroup()


# Match NOA output to each visit

# Find all visits in which an injection might have occurred 
# Exclude visits for glaucoma etc.

# Find just the encounters for the selected eyes from the index date on.
oct_visits <- eye %>% 
  select(PatientID, EyeCode, index_date) %>% 
  inner_join(encounters, by = "PatientID") %>% 
  filter(EncounterDate >= index_date) %>% 
  # Retain only encounters when an OCT was taken
  inner_join(oct_thickness %>% 
               distinct(PatientID, EyeCode, ExamDate)
             , by = c("PatientID", "EyeCode", "EncounterDate" = "ExamDate")) %>% 
  
  # Join the primary outcome, injection or not at each encounter
  left_join(injections %>% 
              select(EyeCode, EncounterID, InjectedDrugDesc),
              by = c("EyeCode", "EncounterID")) %>% 
  # Define the outcome for each visit, aggregating over encounters
  group_by(PatientID, EyeCode, EncounterDate) %>% 
  summarise(injected = any(!is.na(InjectedDrugDesc)), .groups = "drop") %>% 
  
  # Join on the NOA output -  need to know how many were not processed/ineligible here
  left_join(noa, by = c("PatientID", "EyeCode", "EncounterDate" = "OCTDate")) %>% 
  
  # Indicate whether the NOA processing was unsuccessful for a given visit
  # If there were multiple scans taken, all must have been unsuccessful for the visit to be flagged here
  left_join(noa_unsuccessful %>% 
            group_by(PatientID, EyeCode, OCTDate) %>% 
            slice_max(OCTSequence, n=1, with_ties = FALSE) %>% 
            ungroup() %>% 
              select(PatientID, EyeCode, OCTDate, Details), 
              by = c("PatientID", "EyeCode", "EncounterDate" = "OCTDate")) %>% 
  mutate(Details = if_else(!is.na(`Analysis eligibility`), as.character(NA), Details)) %>%
  mutate(across(EncounterDate, as.Date)) %>% 
  
  # Summarise OCT status for each visit
  mutate(OCT_status = case_when(`Analysis eligibility` == 1 ~ "Eligible", 
                                `Analysis eligibility` == 0 ~ "Ineligible",
                                !is.na(Details) ~ "Unsuccessful"))
  
# Just eyes with at least one valid OCT visit
eye_oct <- eye %>% 
  inner_join(oct_visits %>% 
              distinct(PatientID, EyeCode),
             by = c("PatientID", "EyeCode"))



### Sequencing of AVI files ###

# Generate list of AVI files to prioritise
oct_to_process <- oct_visits %>% 
   filter(is.na(`Analysis eligibility`)) %>% 
  select(PatientID, EyeCode, EncounterDate) %>% 
  # Generate the root filenames to search for
  mutate(GenFileName = paste(PatientID, 
                             if_else(EyeCode == "R", "OD", "OS"), 
                             format(EncounterDate, format = "%Y%m%d"),       
                             sep = "_"))

# write_csv(oct_to_process %>% 
#             select(GenFileName),
#           file = find_rstudio_root_file("data", "BIRAX data", "NOA-output", "oct_to_process.csv"))

### Priorities OCTs to process ###
### Execute on NOA machine ###
# 
# # List of OCTs to prioritise
# oct_to_process <- read_csv("F:\\BIRAX\\oct_to_process.csv")
# 
# # List of OCTs available
# oct_available <- list.files(path = "F:\\BIRAX\\HEYEX_Outputs", pattern = ".avi") %>% 
#   enframe(name = NULL, value = "FileName") %>% 
#   mutate(GenFileName = str_remove(FileName, "(_[0-9]{3})*.avi$")) 
#   
# # Find any file NOT in the priority list (anti_join gives all rows from x without a match in y)
# oct_to_exclude <- anti_join(oct_available, oct_to_process, by = "GenFileName") %>% 
#   mutate(FileName = paste0("F:\\BIRAX\\HEYEX_Outputs\\", FileName)) %>% 
#   pull(FileName)
# 
# # Move these to the Excluded scans folder, leaving just those to process in the Heyex folder
# fs::file_move(oct_to_exclude, new_path = "F:\\BIRAX\\Excluded_OCT_Scans")

### End of NOA machine section ###


## Basic descriptive statistics of the cohort

# Patients for whom both eyes were followed-up. 
# i.e. those with fellow eye involvement.

# Visits at which both eyes were imaged
both_eye_visits <- oct_thickness %>% 
  distinct(PatientID, EyeCode, ExamDate) %>% 
  count(PatientID, ExamDate) %>% 
  filter(n == 2) %>% 
  group_by(PatientID) %>% 
  summarise(both_eye_visits = n(),
            oct_start = min(ExamDate),
            oct_end = max(ExamDate)) %>% 
  mutate(oct_both_observed_years = interval(oct_start, oct_end)/dyears())

# Visits at which single eyes were imaged
single_eye_visits <- oct_thickness %>% 
  distinct(PatientID, EyeCode, ExamDate) %>% 
  count(PatientID, ExamDate) %>% 
  filter(n == 1) %>% 
  count(PatientID, name = "single_eye_visits")


# Find all patients with at least one eye with an AMD diagnosis and at least one injection
both_eye_summary <- eye_raw %>% 
  filter(!exclude_age, !exclude_no_AMD, !exclude_RVO, !exclude_DR_DMO, !exclude_no_va, total_injections >=1) %>% 
  distinct(PatientID) %>% 
  as_tibble() %>% 
  
  # Retain just patients where both eyes were imaged simultaneously at least once
  inner_join(both_eye_visits, by = "PatientID") %>% 
  mutate(both_eye_visits = if_else(is.na(both_eye_visits), 0L, both_eye_visits)) %>% 

  # Number of visits where a single eye was imaged
  left_join(single_eye_visits, by = "PatientID") %>% 
  mutate(single_eye_visits = if_else(is.na(single_eye_visits), 0L, single_eye_visits)) %>% 

  # Calculate difference in first injection dates (if both received treatment)
  # To indicate number of events.
  left_join(eye_raw %>% 
              filter(total_injections >= 1) %>% 
              select(PatientID, EyeCode, index_date, total_injections) %>% 
              pivot_wider(values_from = c(index_date, total_injections), names_from = EyeCode) %>% 
              filter(!is.na(index_date_L) | !is.na(index_date_R)) %>% 
              mutate(injection_delay = interval(index_date_L, index_date_R),
                     injection_delay_months = abs(injection_delay/dmonths())), by = "PatientID") %>% 
  # Indicate conversion event in fellow eye
  mutate(fellow_eye_conversion = !is.na(injection_delay_months) & injection_delay_months > 0)
  

