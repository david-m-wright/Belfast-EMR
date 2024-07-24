# Cohort construction for 4T-BIG project 

library(tidyverse)
library(lubridate)

# This is a subset of the main BIRAX cohort
source(rprojroot::find_rstudio_root_file("code/BIRAX-assemble-cohort.R"))

# Prepare dataset

# For each patient calculate the time since the first injection in either eye (first eye index_date)
index_dates_4t <- as.data.table(injection_summary_eye)[, .(index_date = min(index_date)), by = PatientID]

# Dates of clinic visits at the patient level (may be more than one encounter a day so aggregate by date)
visits_4t <- index_dates_4t[
  distinct(encounters_raw, PatientID, EncounterDate), on = c("PatientID"), nomatch = 0]
setkey(visits_4t, PatientID, EncounterDate)
 # Give each patient a visit sequence   
visits_4t[, visit_number := dense_rank(EncounterDate), by = PatientID]  

# Mark the index visit
visits_4t[, index_visit := EncounterDate == index_date]
# Generate a visit sequence number starting at zero for the index visit
visits_4t[, visit_sequence := visit_number - visit_number[index_visit], by = PatientID]
# Calculate number of months since the index visit
visits_4t[, months_since_index := interval(index_date, EncounterDate)/dmonths()]


visits_4t[PatientID == "266211D0-2020-3632-F2B5-EEF012D66060"]

# Exclude patients with neither eye diagnosed with AMD at any point
# Flag eyes with AMD diagnosis at any point
# left_join(AMD_patients %>% 
#             distinct(PatientID) %>% 
#             mutate(exclude_no_AMD = FALSE), by = "PatientID")


# Put all measurements taken in each visit together
visits_4t_summary <- 
visits_4t %>%
  # Indicate if injection occurred.
  left_join(
    dcast(as.data.table(injections_clean)[, .(PatientID, EyeCode = paste0("inj_", EyeCode), EncounterDate = as.POSIXct(EncounterDate), 
                                              injection = 1)],
          PatientID + EncounterDate ~ EyeCode, fill = 0, value.var = "injection"), 
            by = c("PatientID", "EncounterDate")) %>% 

  # Indicate if an AMD diagnosis occurred
  left_join(
    dcast(as.data.table(AMD_patients)[, .(PatientID, EyeCode = paste0("amd_", EyeCode), 
                                          EncounterDate = as.POSIXct(as.Date(DateofDiagnosis_AMD)),
                                          amd = 1)],
          PatientID + EncounterDate ~ EyeCode, fill = 0, value.var = "amd"), 
    by = c("PatientID", "EncounterDate")
  ) %>% 
  
  # Indicate if VA measurement taken
  left_join(
    dcast(va_raw[, .(PatientID, EyeCode = paste0("va_", EyeCode), EncounterDate, va = 1)],
          PatientID + EncounterDate ~ EyeCode, fill = 0, value.var = "va") 
    ,by = c("PatientID", "EncounterDate")) %>% 
  
  # Indicate if OCT volume scan was taken
  left_join(
    dcast(oct_volumes[, .(PatientID, EyeCode = paste0("oct_", EyeCode), EncounterDate = ExamDate, oct = 1)],
          PatientID + EncounterDate ~ EyeCode, fill = 0, value.var = "oct") 
  ,by = c("PatientID", "EncounterDate")) %>% 
  
  # Indicate if OCT multicolour scan was taken
  left_join(
      dcast(oct_multicol[, .(PatientID, EyeCode = paste0("multicol_", EyeCode), EncounterDate = ExamDate, multicol = 1)],
          PatientID + EncounterDate ~ EyeCode, fill = 0, value.var = "multicol") 
  ,by = c("PatientID", "EncounterDate")) %>% 
  
  left_join(
    # Indicate if fluorescein angiogram was taken
    dcast(oct_fa[, .(PatientID, EyeCode = paste0("fa_", EyeCode), EncounterDate = ExamDate, fa = 1)],
          PatientID + EncounterDate ~ EyeCode, fill = 0, value.var = "fa") 
  ,by = c("PatientID", "EncounterDate")) %>% 
  
  mutate(across(matches("_L$|_R$"), ~if_else(is.na(.), 0, .)))



# visits_4t_summary[index_visit == TRUE & inj_L == 1] %>% 
#   count(amd_L)
#   count(va_L, oct_L, multicol_L, fa_L)
#   # count(inj_L, inj_R)
#   count(inj_L, inj_R, fa_L, fa_R)

# Active nAMD at baseline
# This will always contain the first eye or both if a bilateral case
namd_baseline <-  visits_4t_summary[index_visit == TRUE & 
                                        ((inj_R == 1 & amd_R == 1)|
                                           (inj_L == 1 & amd_L == 1))][,
                                      # Indicate eyes that are just missing multicolour images at baseline
                                    c("no_multicol_L",
                                      "no_multicol_R") := list((va_L ==1 & oct_L == 1 & fa_L==1),
                                                                 (va_R == 1 & oct_R == 1 & fa_R==1))
                                    ]
namd_baseline[inj_L == TRUE,] %>% 
  count(va_L, oct_L, multicol_L, fa_L, no_multicol_L) %>% 
  adorn_totals("row")

namd_baseline[inj_R == TRUE,] %>% 
  count(va_R, oct_R, multicol_R, fa_R, no_multicol_R) %>% 
  adorn_totals("row")

# Imaging by months since baseline
# Available imaging after >5 years followup
namd_baseline[inj_L == TRUE & no_multicol_L == TRUE,.(PatientID)][
  visits_4t_summary, on = "PatientID", nomatch =0][
    months_since_index >=60][,  
    # months_since_index := floor(months_since_index)][,
      .(va_L = sum(va_L)>0,
          oct_L = sum(oct_L)>0,
          multicol_L = sum(multicol_L)>0,
          fa_L = sum(fa_L)>0), 
      by = PatientID
      # by = .(PatientID, months_since_index)
    ] %>% count(va_L, oct_L, multicol_L, fa_L) %>% 
  adorn_totals("row")


namd_baseline[inj_R == TRUE & no_multicol_R == TRUE,.(PatientID)][
  visits_4t_summary, on = "PatientID", nomatch =0][
    months_since_index >=60][,  
                             # months_since_index := floor(months_since_index)][,
                             .(va_R = sum(va_R)>0,
                               oct_R = sum(oct_R)>0,
                               multicol_R = sum(multicol_R)>0,
                               fa_R = sum(fa_R)>0), 
                             by = PatientID
                             # by = .(PatientID, months_since_index)
    ] %>% count(va_R, oct_R, multicol_R, fa_R) %>% 
  adorn_totals("row")



# Available imaging after 7 years followup
namd_baseline[inj_L == TRUE & no_multicol_L == TRUE,.(PatientID)][
  visits_4t_summary, on = "PatientID", nomatch =0][
    months_since_index >=84][,  
                             # months_since_index := floor(months_since_index)][,
                             .(va_L = sum(va_L)>0,
                               oct_L = sum(oct_L)>0,
                               multicol_L = sum(multicol_L)>0,
                               fa_L = sum(fa_L)>0), 
                             by = PatientID
                             # by = .(PatientID, months_since_index)
    ] %>% count(va_L, oct_L, multicol_L, fa_L) %>% 
  adorn_totals("row")


namd_baseline[inj_R == TRUE & no_multicol_R == TRUE,.(PatientID)][
  visits_4t_summary, on = "PatientID", nomatch =0][
    months_since_index >=84][,  
                             # months_since_index := floor(months_since_index)][,
                             .(va_R = sum(va_R)>0,
                               oct_R = sum(oct_R)>0,
                               multicol_R = sum(multicol_R)>0,
                               fa_R = sum(fa_R)>0), 
                             by = PatientID
                             # by = .(PatientID, months_since_index)
    ] %>% count(va_R, oct_R, multicol_R, fa_R) %>% 
  adorn_totals("row")


# Count visits and images for patients with followup >7 years after index date
patients_7yr <- distinct(namd_baseline[(inj_R == T & oct_R == 1)|(inj_L == T & oct_L == 1),.(PatientID)][
  visits_4t_summary, on = "PatientID", nomatch =0][
  months_since_index >=84 & (oct_L == 1 | oct_R == 1), .(PatientID)])[
    visits_4t_summary, on = "PatientID", nomatch = 0
  ][, months_since_index := floor(months_since_index)][,
    .(visits = .N,
      injections = sum(inj_L+inj_R),
      va_measurements = sum(va_L + va_R),
      oct = sum(oct_L + oct_R),
      multicolour = sum(multicol_L + multicol_R),
      fa = sum(fa_L + fa_R)),
      by = .(months_since_index, index_visit)][
        order(months_since_index)
      ] 

fwrite(patients_7yr, find_rstudio_root_file("data", "output", "7yr_followup.csv"))







