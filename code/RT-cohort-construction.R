# Cohort construction for Rosetrees Trust project on prediction fellow eye conversion to AMD

library(tidyverse)
library(lubridate)

source(rprojroot::find_rstudio_root_file("code/BIRAX-assemble-cohort.R"))


# Prepare dataset

# For each visual acuity measurement, calculate the time since the first injection in the first eye (first eye index_date)
# Differs to BIRAX which had eye specific index dates
va_raw_rt <- as.data.table(injection_summary_eye)[eye_order == 1, .(PatientID, index_date, EyeCode_1 = EyeCode)][ 
  # Join injection_summary_eye and visual_acuity tables
  as.data.table(visual_acuity), .(PatientID, EyeCode, EyeCode_1, index_date, EncounterDate, va_logmar, va_etdrs, va_category_snellen), on = .(PatientID)][
    # Calculate months since index date
    , months_since_index:=as.interval(index_date, EncounterDate)/dmonths()][
      # Calculate years ince index date
      , years_since_index := as.interval(index_date, EncounterDate)/dyears()][  
        # Mark baseline measurements (closest measurement to baseline)  
        , baseline := abs(months_since_index) == min(abs(months_since_index)), by = .(PatientID, EyeCode)][
           # Add the eye order - eyes with a tied first injection have eye order set to NA
          , eye_order := case_when(EyeCode == EyeCode_1 & !is.na(EyeCode_1) ~ 1,
                                  EyeCode != EyeCode_1 & !is.na(EyeCode_1) ~ 2,)]



# Details of all the examinations that were extracted
# Note that there were multiple ExamIDs and SeriesIDs (scans within exams) on some dates
oct_details <- fread(file = file.path(file_path, "OCT_ImageVariables.txt"))

# List of all dates on which an OCT volume scan was taken 
# When there are multiple exams on the same day, select the last (matches criteria used for raw NOA data)
# Assumes there was a problem with the earlier scan
oct_volumes <- unique(oct_details[ModalityType == "Volume" & ModalityProcedure == "IR_OCT" & ImageType == "OCT",], by = c("PatientID", "EyeCode", "ExamDate"), fromLast = TRUE)


## Generate a history of OCT scans for each eye 

# Find date of first anti-VEGF injection for each eye (index dates)
index_dates_rt <- injection_summary_eye %>% 
  filter(eye_order != 1.5) %>% 
  select(PatientID, EyeCode, index_date, eye_order) %>% 
  pivot_wider(names_from = eye_order, values_from = c("EyeCode", "index_date")) %>% 
  as.data.table()

# Add a variable to each data.table to make joins easier to navigate (otherwise one of the join columns will disappear)
# This is especially important for rolling joins
# Date format for the join columns must be identical
oct_volumes[, visit_date := as.Date(ExamDate)]
index_dates_rt[, join_date := as.Date(index_date_1)]

# Set the join keys and sort the data - this must be done before determining OCT sequences
setkey(oct_volumes, PatientID, visit_date)
setkey(index_dates_rt, PatientID, join_date)


# Mark each OCT scan with a visit number, with one being the baseline visit on the day of the first injection (index_date)
oct_history_raw_rt <-
  index_dates_rt[, .(PatientID, index_date_1)][oct_volumes, on = "PatientID"][, c("visit_number", "months_since_index") := list(
    # Visit number - if an OCT was missed for one eye that visit will be missed from that eye sequence
    dense_rank(visit_date),
    # Time since first injection - negative if the visit was before the index date - so visit 1 may be before index_visit for the first eye
    interval(index_date_1, visit_date) /
      dmonths()),
    by = "PatientID"][,
                      # Mark the closest OCT visit to the index date up to 2 weeks prior to the index date
                      # index_visit := .I == .I[which.min(abs(months_since_index))],    by = "PatientID"][
                        index_visit := abs(months_since_index) == min(abs(months_since_index)),    by = "PatientID"][
                        # OCT visit sequence starting at zero at the index visit. Negative values were OCTs taken before the index_visit
                        , visit_sequence := visit_number - visit_number[index_visit], by = "PatientID"]
setkey(oct_history_raw_rt, PatientID, visit_date)


## Patient level information

patients_raw_rt  <- patients %>% 
    select(SiteID, PatientID, PerturbedDateofBirth, PerturbedCurrentAge, PerturbedAgeatDeath, Gender) %>% 
    
    # Exclude patients with neither eye diagnosed with AMD at any point
    left_join(AMD_patients %>% 
                distinct(PatientID) %>% 
                mutate(exclude_no_AMD = FALSE), by = "PatientID") %>% 
    mutate(exclude_no_AMD = if_else(is.na(exclude_no_AMD), TRUE, exclude_no_AMD)) %>% 
  

    # Exclude patients with neither eye injected with anti-VEGF at any point
    left_join(injection_summary_eye %>% 
            distinct(PatientID) %>% 
            mutate(exclude_no_injections = FALSE), by = "PatientID") %>% 
    mutate(exclude_no_injections = if_else(is.na(exclude_no_injections), TRUE, exclude_no_injections)) %>% 

    
    # Exclude patients with both eyes converting on the same date
    left_join(injection_summary_eye %>%
        filter(eye_order == 1.5) %>%
        distinct(PatientID) %>%
        mutate(exclude_bilateral_index = TRUE), by = "PatientID") %>%
    mutate(exclude_bilateral_index = if_else(is.na(exclude_bilateral_index), FALSE, exclude_bilateral_index)) %>%


    # Find date of first anti-VEGF injection for each eye (index dates)
    left_join(index_dates_rt,  by = "PatientID")  %>% 
    mutate(EyeCode_2 = if_else(EyeCode_1 == "L", "R", "L")) %>% 
    
    
    # Join VA measured on the person level index date or up to 2 weeks prior as the baseline
    left_join(va_raw_rt %>%
              filter(baseline, months_since_index > -0.5 & months_since_index <=0) %>% 
                select(PatientID, eye_order, va_logmar, va_etdrs) %>% 
                pivot_wider(names_from = eye_order, values_from = c("va_logmar", "va_etdrs")),
            by = "PatientID") %>% 
    mutate(exclude_no_va = is.na(va_logmar_1) | is.na(va_logmar_2) | is.na(va_etdrs_1) | is.na(va_etdrs_2)) %>% 
  
    # Exclude patients with neither eye injected with anti-VEGF at any point
    mutate(exclude_no_injections = is.na(index_date_1) & is.na(index_date_2)) %>% 
  
    # Exclude those with RVO or DMO/DR
    left_join(DMO_DR_patients %>% 
            distinct(PatientID) %>% 
            mutate(exclude_DR_DMO = TRUE), by = "PatientID") %>% 
    mutate(exclude_DR_DMO = if_else(is.na(exclude_DR_DMO), FALSE, exclude_DR_DMO)) %>% 
      left_join(RVO_patients %>% 
                  distinct(PatientID) %>% 
                  mutate(exclude_RVO = TRUE), by = "PatientID") %>% 
      mutate(exclude_RVO = if_else(is.na(exclude_RVO), FALSE, exclude_RVO)) %>% 


  # Exclude those without OCT of both eyes on the index date or within 2 weeks prior to the index
  left_join(oct_history_raw_rt[index_visit & months_since_index >-0.5 & months_since_index <= 0, .N, by = "PatientID"][, exclude_no_oct := N!=2], by = "PatientID") %>% 
  mutate(exclude_no_oct = if_else(is.na(exclude_no_oct), TRUE, exclude_no_oct)) %>% 

  # Put in the OCT sequence start (start of follow up) and end dates for patient
  left_join(oct_history_raw_rt[(index_visit & months_since_index > -0.5) | months_since_index >= 0, .(fup_start = min(visit_date), oct_end = max(visit_date), oct_visits = max(visit_sequence)+1), by = "PatientID"], by = "PatientID") %>% 
  
  # Exclude those with only a single OCT visit (not possible to assess outcome as not being regularly monitored)
  mutate(exclude_single_oct = oct_visits == 1) %>% 
  
    # Add the end of follow-up and outcome
  mutate(fup_end = coalesce(index_date_2, oct_end),
         conversion = !is.na(index_date_2),
         fup_months = interval(fup_start, fup_end)/dmonths()) %>% 
  
  # Calculate age at start of follow-up
  mutate(index_age =  interval(PerturbedDateofBirth, fup_start)/dyears()) %>% 
  
  # Exclude those aged <50 or younger on date of first injection (index date) (or no injections)
  mutate(exclude_age = index_age < 50 & !is.na(index_age))
  

# Check exclusions
  patients_raw_rt %>% 
      count(across(matches("exclude"))) %>% 
      arrange(n)
    

# Apply the exclusion criteria
patients_rt <- patients_raw_rt %>% 
  filter(!exclude_no_AMD,
         !exclude_no_injections,
         !exclude_bilateral_index,
         !exclude_no_va,
         !exclude_age,
         !exclude_DR_DMO,
         !exclude_RVO,
         !exclude_no_oct,
         !exclude_single_oct) %>% 
  select(-N, -matches("exclude"))



# The OCT history just for cohort members 
oct_history_rt <- patients_rt[, .(PatientID, fup_start, fup_end, conversion, EyeCode_1)][oct_history_raw_rt, on = c(PatientID = "PatientID"), nomatch = 0][
  visit_date >= fup_start & visit_date <= fup_end,][
    , eye_order := if_else(EyeCode == EyeCode_1, 1, 2)
  ]

