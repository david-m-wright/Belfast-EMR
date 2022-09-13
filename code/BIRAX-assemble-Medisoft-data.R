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
         # Note that the 2DP logmar is rounded to 1DP before conversion to ETDRS
         logmar_to_etdrs = as.character(coalesce(va(Logmar2DBestMeasure, from = "logmar", to = "etdrs"),
                                                 va(Logmar1DBestMeasure, from = "logmar", to = "etdrs"))),
         va_etdrs = as.numeric(coalesce(cleaned_etdrs, snellen_to_etdrs, logmar_to_etdrs))) %>% 
  
  mutate(cleaned_logmar = as.character(coalesce(va(Logmar2DBestMeasure, from = "logmar", to = "logmar"), 
                                                va(Logmar1DBestMeasure, from = "logmar", to = "logmar"))),
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
                                             na_level = "Missing")) %>% 
  # Take best VA measurement on a given day
  group_by(PatientID, EyeCode, EncounterDate) %>% 
  slice_min(va_logmar, n=1, with_ties = FALSE) %>% 
  ungroup()


# Injections
injections_raw <- fread(file.path(file_path, "BELInjections.txt")) %>% 
  # Exclude non anti-VEGF injections
  mutate(exclude_not_antiVEGF = AntiVEGFInjection != 1) %>% 
  # Exclude a small number of duplicate encounters on the same date so only one injection per eye per date
  group_by(PatientID, EyeCode, EncounterDate, exclude_not_antiVEGF) %>% 
  slice_head(n=1) %>% 
  ungroup()
  
injections_clean <- injections_raw %>% 
  filter(!exclude_not_antiVEGF) %>% 
  select(PatientID, EyeCode, EncounterDate, InjectedDrugDesc, EncounterID) #%>% 
  # # Add visual acuity at the selected encounter
  # left_join(visual_acuity %>%
  #             select(PatientID, EyeCode, EncounterDate, va_etdrs, va_logmar),
  #           by = c("PatientID", "EyeCode", "EncounterDate")) 

# OCT thickness maps
oct_thickness <- fread(file.path(file_path, "OCT_ThicknessMetrics.txt")) %>% 
  # Select a single scan for each eye and date
  group_by(PatientID, EyeCode, ExamDate) %>% 
  # ScanSeq is the latest exam for that eye on that date
  slice_min(order_by = ScanSeq, with_ties = FALSE) %>% 
  ungroup()


