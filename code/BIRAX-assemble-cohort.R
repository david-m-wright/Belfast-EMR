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




### Prepare time series of measurements over time ###

# These analyses are all based on the injection history for each eye, starting at the date of the first injection

# How do we know these eyes are treatment naive? Especially those at the start of the dataset? *********
# Possible sensitivity analysis excluding treatment start in the first 6 months of the dataset?

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
  .[as.data.table(visual_acuity), .(PatientID, EyeCode, index_date, EncounterDate, va_logmar, va_etdrs, va_category_snellen), on = .(PatientID, EyeCode)] %>% 
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
# Baseline dates must match exactly
snapshots <- data.table(months_since_index = c(0, 12, 24, 36, 48, 60, 84, 120), 
                        follow_up = c(0, 12, 24, 36, 48, 60, 84, 120), 
                        tolerance = c(0, 2, 2, 2, 2, 2,2,2),
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

# Retain only those with a VA measurement close enough to the index date to be called a baseline (now set to exact match on baseline)
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


## OCT thickness metrics


# For each OCT thickness measurement, calculate the time since the first injection (index_date)
oct_thickness_raw <- injection_summary_eye %>% 
  as.data.table() %>% 
  # Join injection_summary_eye and oct_thickness tables
  .[as.data.table(oct_thickness), on = .(PatientID, EyeCode)] %>% 
  # Retain OCT measurements from the index date onward
  .[ExamDate >= index_date] %>% 
  .[, months_since_index:=as.interval(index_date, ExamDate)/dmonths()] 


# Calculate changes in thickness measurements from baseline 
# Could do this with a set of baseline columns and a set of change columns
# ss <- oct_thickness_raw %>% 
#   .[, .(PatientID, EyeCode, ExamDate, OuterSuperior_Volume, OuterInferior_Volume)] 
#   cnames <- c("OuterSuperior_Volume", "OuterInferior_Volume")
#   ss[, .SDcols = cnames,  (paste0(cnames, "_change")) := lapply(.SD, function(x) x - x[1]), by=list(PatientID, EyeCode)] %>% 
#     print()


## NOA fluid measurements
  # For each OCT thickness measurement, calculate the time since the first injection (index_date)
  fluid_raw <- injection_summary_eye %>% 
    as.data.table() %>% 
    # Join injection_summary_eye and oct_thickness tables
    .[as.data.table(noa), on = .(PatientID, EyeCode)] %>% 
    # Retain OCT measurements from the index date onward
    .[OCTDate >= index_date] % >% 
    .[, months_since_index:=as.interval(index_date, OCTDate)/dmonths()] 
  
  
  

# Tabulate the measures taken at each visit


# Single VA measurement for each date - but need to choose the one that best matches the injection EncounterID???

# PatientID | EyeCode | Date | VA | Injection | OCT thickness | Fluid

# Round observations to months if dates do not exactly match?
observations <- 
  # visual_acuity %>%
  # transmute(PatientID, EyeCode, EncounterDate, va = TRUE) %>%
  
  va_history_clean %>% 
  filter(snapshot == TRUE & follow_up == 0) %>% 
  transmute(PatientID, EyeCode, EncounterDate, va = TRUE) %>%
  
  full_join(
    injections_clean %>%
      transmute(PatientID, EyeCode, EncounterDate, inj = TRUE),
    by = c("PatientID", "EyeCode", "EncounterDate")
  ) %>%
  
  full_join(
    # Note that there are a lot of missing values in the thickness measurements
    oct_thickness_raw %>%
      transmute(PatientID, EyeCode, ExamDate, oct = TRUE),
    by = c("PatientID", "EyeCode", "EncounterDate" = "ExamDate")
  ) %>%
  
  full_join(
    fluid_raw %>%
      transmute(PatientID, EyeCode, OCTDate, fluid = TRUE),
    by = c("PatientID", "EyeCode", "EncounterDate" = "OCTDate")
  )


# Timing of observations does not match well
observations %>% count(va, inj, oct, fluid)


# Use rolling joins in data.table to roll values forward
# Approach for descriptive statistics - roll forward or report with missingness?
# Report with missingness

# Missing values for the covariates?


# Assemble eye level dataset (single row per eye) 

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
  
  # Add VA measured on the index date
  left_join(va_history_clean %>%
              transmute(PatientID, EyeCode, EncounterDate, va_logmar, "ETDRS letters" = va_etdrs),
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
    exclude_no_intervals = total_intervals == 0,
    
    # Those without a baseline VA measurement
    exclude_no_va = is.na(va_logmar)) 



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
oct_history <- oct_thickness_raw %>%  
  inner_join(eye %>% 
               select(PatientID, EyeCode, index_date), by = c("PatientID", "EyeCode", "index_date")) %>% 
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



# Note that the scan sequences from the OCT thickness table are ordered differently from the NOA filenames table. 
# Therefore it may be difficult to directly match the OCT scan used to generate each metric.



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


# 
# # Are there fellow eye images?
# oct_thickness %>% 
#   distinct(PatientID, ExamDate, EyeCode) %>% 
#   count(PatientID, ExamDate) %>% 
#   count(n)
# # Yes - approx 85% of dates
# 
# # Are there fellow eye fluid measurements?
# noa %>% 
#   distinct(PatientID, OCTDate, EyeCode) %>% 
#   count(PatientID, OCTDate) %>% 
#   count(n)
# 
