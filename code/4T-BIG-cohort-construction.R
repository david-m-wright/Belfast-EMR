# Cohort construction for 4T-BIG and fibrosis projects

library(tidyverse)
library(lubridate)
library(fs)

# This is a subset of the main BIRAX cohort
source(rprojroot::find_rstudio_root_file("code/BIRAX-assemble-cohort.R"))

## Prepare dataset

# Details of all the examinations that were extracted with imaging
# Note that there were multiple ExamIDs and SeriesIDs (scans within exams) on some dates
oct_details <- fread(file = file.path(file_path, "OCT_ImageVariables.txt"))
oct_details[, c("StartCoordX", "StartCoordY") := lapply(.SD, as.numeric), .SDcols = c("StartCoordX", "StartCoordY")]
setkey(oct_details, FilePath)

# List of all dates on which an OCT volume scan was taken 
# When there are multiple exams on the same day, select the last (matches criteria used for raw NOA data)
# Assumes there was a problem with the earlier scan
oct_volumes <- unique(oct_details[ModalityType == "Volume" & ModalityProcedure == "IR_OCT" & ImageType == "OCT",], by = c("PatientID", "EyeCode", "ExamDate"), fromLast = TRUE)[, 
                # Add a volume identifier
                volume_id := str_replace(FilePath, "(.*\\\\Volume\\\\.*)(-.*)", "\\1")]

# List of all dates where an OCT multicolour photograph was taken (again, take last scan if there were multiple on the same date)
oct_multicol <- unique(oct_details[ModalityType == "Single" & ModalityProcedure == "MC" & ImageType == "ANGIO",], by = c("PatientID", "EyeCode", "ExamDate"), fromLast = TRUE)

# List of all dates where a fluorescein angiogram (en face) was taken (again, take last scan if there were multiple on the same date)
oct_fa <- unique(oct_details[(ModalityProcedure == "FA" & ModalityType == "Single"),], by = c("PatientID", "EyeCode", "ExamDate"), fromLast = TRUE)

# For each patient calculate the time since the first injection in either eye (first eye index_date)
index_dates_4t <- as.data.table(injection_summary_eye)[, .(index_date = min(index_date)), by = PatientID]

# Dates of clinic visits at the patient level (may be more than one encounter a day so aggregate by date)
visits_4t <- index_dates_4t[
  distinct(encounters_raw, PatientID, EncounterDate), on = c("PatientID"), nomatch = 0]
setkey(visits_4t, PatientID, EncounterDate)
 # Give each patient a visit ID   
visits_4t[, visit_number := dense_rank(EncounterDate), by = PatientID]  

# Mark the index visit
visits_4t[, index_visit := EncounterDate == index_date]
# Generate a visit sequence number starting at zero for the index visit
visits_4t[, visit_sequence := visit_number - visit_number[index_visit], by = PatientID]
# Calculate number of months since the index visit
visits_4t[, months_since_index := interval(index_date, EncounterDate)/dmonths()]

# Note that some visits did not have any imaging
# encounters_raw[PatientID == "FC9B5F4A-EE14-2E67-8E4E-EB5C162C6D71" & EncounterDate >= "2017-07-22"][order(EncounterDate)]


# Put all measurements taken in each visit together
visits_4t_summary <- visits_4t %>%
  
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
    # Indicate if fluorescein angiogram was taken (either OCT or 2D)
    dcast(oct_fa[, .(PatientID, EyeCode = paste0("fa_", EyeCode), EncounterDate = ExamDate, fa = 1)],
          PatientID + EncounterDate ~ EyeCode, fill = 0, value.var = "fa") 
  ,by = c("PatientID", "EncounterDate")) %>% 
  
  mutate(across(matches("_L$|_R$"), ~if_else(is.na(.), 0, .))) %>% 

  # Drop visits where not VA, injection or imaging took place (may have been booking visits)
  # This will break the visit sequences for each patient
  filter(!(inj_L == 0 & inj_R == 0 & 
             va_L == 0 & va_R == 0 & 
             oct_L == 0 & oct_R == 0 & 
             multicol_L == 0 & multicol_R == 0 & 
             fa_L == 0 & fa_R == 0))
 


# Inclusion criterion: active nAMD at baseline
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

# Distribution of imaging by months since baseline
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



### Cohorts ####

## 4T-BIG study

# Only include patients with a visit at exactly 3 years

patients_3yr_exact <- distinct(namd_baseline[(inj_R == T & oct_R == 1)|(inj_L == T & oct_L == 1),.(PatientID)][
  visits_4t_summary, on = "PatientID", nomatch =0][
    months_since_index >=36 & months_since_index < 37 & (oct_L == 1 | oct_R == 1), .(PatientID)])[
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
                                                         ][
                                                           months_since_index <=36 ] 

# fwrite(patients_3yr_exact, find_rstudio_root_file("data", "output", "3yr-followup-exact.csv"))


## Fibrosis project

# List of patient visits with imaging status 
# Only include patients with a visit at exactly 6 years
patients_6yr_exact_list <- distinct(namd_baseline[(inj_R == T & oct_R == 1)|(inj_L == T & oct_L == 1),.(PatientID)][
  # Joined to visit summary
  visits_4t_summary, on = "PatientID", nomatch =0][
    # Just patients with an OCT visit in the 72nd month (double)
    months_since_index >=72 & months_since_index < 73 & (oct_L == 1 | oct_R == 1), .(PatientID)])[
      visits_4t_summary, on = "PatientID", nomatch = 0
    ]
# fwrite(patients_6yr_exact_list, "F:\\Fibrosis\\visit-list.csv"))


# Patients with a valid baseline (injection and OCT in either eye)
patients_6yr_exact <- patients_6yr_exact_list[
      # Convert time since index to integer
      , months_since_index := floor(months_since_index)][
        # Calculate number of visits/images by month
        ,
        .(visits = .N,
          injections = sum(inj_L+inj_R),
          va_measurements = sum(va_L + va_R),
          oct = sum(oct_L + oct_R),
          multicolour = sum(multicol_L + multicol_R),
          fa = sum(fa_L + fa_R)),
        by = .(months_since_index, index_visit)][
          order(months_since_index)
        ][
          # Retain only first 7 years of followup
          months_since_index <=72 ] 
# fwrite(patients_6yr_exact, find_rstudio_root_file("data", "output", "6yr-followup-exact.csv"))


# 
# distinct(namd_baseline[(inj_R == T & oct_R == 1)|(inj_L == T & oct_L == 1),.(PatientID)][
#        visits_4t_summary, on = "PatientID", nomatch =0][
#              months_since_index >=84 & months_since_index < 85 & (oct_L == 1 | oct_R == 1), .(PatientID)])[
#                    visits_4t_summary, on = "PatientID", nomatch = 0
#                ][PatientID == "FC9B5F4A-EE14-2E67-8E4E-EB5C162C6D71"]
# 
# encounters_raw[PatientID == "FC9B5F4A-EE14-2E67-8E4E-EB5C162C6D71" & EncounterDate >= "2017-07-22"][order(EncounterDate)]


# 
# 
# distinct(namd_baseline[(inj_R == T & oct_R == 1)|(inj_L == T & oct_L == 1),.(PatientID)][
#   visits_4t_summary, on = "PatientID", nomatch =0][
#     months_since_index >=84 & months_since_index < 85 & (oct_L == 1 | oct_R == 1), .(PatientID)])[
#       visits_4t_summary, on = "PatientID", nomatch = 0
#     ][
#       months_since_index >=84 & months_since_index <85 ] 
# 
# encounters_raw[PatientID == "0AE95C3D-873D-B865-346F-2502DFCFD372" & EncounterDate >= "2018-06-02" & EncounterDate <= "2018-06-16",]
# 



### Extract full imaging dataset for the selected cohort ###

# All injections as data.table for meta data
injections_dt <- data.table(injections_clean)[, .(PatientID, EyeCode, ExamDate = EncounterDate, injection = 1)]

# OCT scan metadata
oct_volume_metadata <- injections_dt[
  oct_volumes[patients_6yr_exact_list[, .(PatientID, EncounterDate, visit_sequence, months_since_index)], on = c("PatientID","ExamDate" = "EncounterDate"), nomatch = 0],
  on = c("PatientID", "EyeCode", "ExamDate")][, injection := if_else(is.na(injection), 0, injection)]

# Find all the BMP files associated with the selected OCT volumes
# Get the names and paths to the OCT volumes
fnames <- oct_volume_metadata %>% 
  select(FilePath) %>% 
  # slice(1:5) %>% 
  separate_wider_regex(cols = "FilePath", patterns = c(patientpath = ".*", volpath =  "\\\\Volume\\\\", volname = ".*", "-", suffix = ".*"), cols_remove = FALSE)

# Find all the slice filenames by searching on the volume identifier (ignoring the slice specific identifiers) (6mins)
oct_slices_fnames <- map2(.x = paste0(fnames$patientpath, fnames$volpath), .y = fnames$volname, ~list.files(path = paste0("F:\\BIRAX_ProcessedOutputs", .x), pattern = .y, full.names = TRUE), .progress = TRUE) %>% 
  list_c()
oct_slices <- data.table(
   FilePath = str_replace(oct_slices_fnames, "F:\\\\BIRAX_ProcessedOutputs", "")
  # Find volume filename for matching back to volume details
 )[, volume_id := str_replace(FilePath, "(.*\\\\Volume\\\\.*)(-.*)", "\\1")]
setkey(oct_slices, FilePath, volume_id)


out_dir <- find_rstudio_root_file("data", "output", "Fibrosis")
# out_dir <- "E:\\Fibrosis"
# 
# 
# ## Copy imaging to the target directory, maintaining the source directory structure
# 
# # Create a directory for each individual
# dir_create(paste0(out_dir, "\\Imaging", unique(patients_6yr_exact_list$PatientID)))
# 
# # Within each patient, create a directory for the OCT slices
# dir_create(paste0(out_dir, "\\Imaging", unique(fnames$patientpath), "\\Volume"))
# 
# # Where present, create a directory for the multicolour enface images
oct_multicol_metadata <- oct_multicol[patients_6yr_exact_list[, .(PatientID, EncounterDate, visit_sequence, months_since_index)], on = c("PatientID","ExamDate" = "EncounterDate"), nomatch = 0]

oct_multicol_metadata[]
as.data.table(injections_clean)[, .(PatientID, EyeCode, EncounterDate)]

# [, injected := 1]

# dir_create(paste0(out_dir, "\\Imaging", unique(str_extract(oct_multicol_metadata$FilePath,  ".*\\\\Single"))))
# 
# # Where present (and not already created for multicolour), create a directory for FA enface images 
# # (All are numbered slice 1)
oct_fa_metadata <- oct_fa[patients_6yr_exact_list[, .(PatientID, EncounterDate, visit_sequence, months_since_index)], on = c("PatientID","ExamDate" = "EncounterDate"), nomatch = 0]



# dir_create(paste0(out_dir, "\\Imaging", unique(str_extract(oct_fa_metadata$FilePath,  ".*\\\\Single"))))
# 
# 
# # Populate with filtered scans
# # OCT
# file_copy(oct_slices$filename, paste0(out_dir, "\\Imaging", oct_slices$FilePath))
# # Multicolour (8 mins - 12Gb)
# file_copy(paste0("F:\\BIRAX_ProcessedOutputs", oct_multicol_metadata$FilePath), paste0(out_dir, "\\Imaging", oct_multicol_metadata$FilePath))
# # FA
# file_copy(paste0("F:\\BIRAX_ProcessedOutputs", oct_fa_metadata$FilePath), paste0(out_dir, "\\Imaging", oct_fa_metadata$FilePath))
# 
# 
# 
# ## Output the metadata 
# 
# dir_create(paste0(out_dir, "\\Metadata"))
# # For OCT, output a row for each slice.
oct_details[oct_slices[oct_volume_metadata[, .(volume_id, visit_sequence, months_since_index, injection)], on = "volume_id"], on = "FilePath", nomatch = 0]# %>% 
   fwrite(paste0(out_dir, "\\Metadata", "\\OCT-slice-metadata.csv"))
# # # For multiocolour and FA, a single row for each enface image
injections_dt[oct_multicol_metadata, on = c("PatientID", "EyeCode", "ExamDate")][, injection := if_else(is.na(injection), 0, injection)] #%>% 
  fwrite(paste0(out_dir, "\\Metadata", "\\multicol-metadata.csv"))
injections_dt[oct_fa_metadata, on = c("PatientID", "EyeCode", "ExamDate")][, injection := if_else(is.na(injection), 0, injection)] #%>% 
 fwrite(paste0(out_dir, "\\Metadata", "\\fa-metadata.csv"))
# # 
# # 14:53 multicol
# # 15:11 OCT
# 






### Extract clinical data for the selected cohort ###




####################################

## Other potential cohorts

# Only include patients with a visit at exactly 4 years
# Patietns with a valid baseline (injection and OCT in either eye)
patients_4yr_exact <- distinct(namd_baseline[(inj_R == T & oct_R == 1)|(inj_L == T & oct_L == 1),.(PatientID)][
  # Joined to visit summary
  visits_4t_summary, on = "PatientID", nomatch =0][
    # Just patients with a visit in the 48th month (double)
    months_since_index >=48 & months_since_index < 49 & (oct_L == 1 | oct_R == 1), .(PatientID)])[
      visits_4t_summary, on = "PatientID", nomatch = 0
    ][
      # Convert time since index to integer
      , months_since_index := floor(months_since_index)][
        # Calculate number of visits/images by month
        ,
        .(visits = .N,
          injections = sum(inj_L+inj_R),
          va_measurements = sum(va_L + va_R),
          oct = sum(oct_L + oct_R),
          multicolour = sum(multicol_L + multicol_R),
          fa = sum(fa_L + fa_R)),
        by = .(months_since_index, index_visit)][
          order(months_since_index)
        ][
          # Retain only first 7 years of followup
          months_since_index <=48 ] 
# fwrite(patients_4yr_exact, find_rstudio_root_file("data", "output", "4yr-followup-exact.csv"))


# Only include patients with a visit at exactly 5 years
# Patietns with a valid baseline (injection and OCT in either eye)
patients_5yr_exact <- distinct(namd_baseline[(inj_R == T & oct_R == 1)|(inj_L == T & oct_L == 1),.(PatientID)][
  # Joined to visit summary
  visits_4t_summary, on = "PatientID", nomatch =0][
    # Just patients with a visit in the 60th month (double)
    months_since_index >=60 & months_since_index < 61 & (oct_L == 1 | oct_R == 1), .(PatientID)])[
      visits_4t_summary, on = "PatientID", nomatch = 0
    ][
      # Convert time since index to integer
      , months_since_index := floor(months_since_index)][
        # Calculate number of visits/images by month
        ,
        .(visits = .N,
          injections = sum(inj_L+inj_R),
          va_measurements = sum(va_L + va_R),
          oct = sum(oct_L + oct_R),
          multicolour = sum(multicol_L + multicol_R),
          fa = sum(fa_L + fa_R)),
        by = .(months_since_index, index_visit)][
          order(months_since_index)
        ][
          # Retain only first 7 years of followup
          months_since_index <=60 ] 
# fwrite(patients_5yr_exact, find_rstudio_root_file("data", "output", "5yr-followup-exact.csv"))


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

# fwrite(patients_7yr, find_rstudio_root_file("data", "output", "7yr_followup.csv"))


# Only include patients with a visit at exactly 7 years
# Patietns with a valid baseline (injection and OCT in either eye)
patients_7yr_exact <- distinct(namd_baseline[(inj_R == T & oct_R == 1)|(inj_L == T & oct_L == 1),.(PatientID)][
  # Joined to visit summary
  visits_4t_summary, on = "PatientID", nomatch =0][
    # Just patients with a visit in the 84th month (double)
    months_since_index >=84 & months_since_index < 85 & (oct_L == 1 | oct_R == 1), .(PatientID)])[
      visits_4t_summary, on = "PatientID", nomatch = 0
    ][
      # Convert time since index to integer
      , months_since_index := floor(months_since_index)][
        # Calculate number of visits/images by month
        ,
        .(visits = .N,
          injections = sum(inj_L+inj_R),
          va_measurements = sum(va_L + va_R),
          oct = sum(oct_L + oct_R),
          multicolour = sum(multicol_L + multicol_R),
          fa = sum(fa_L + fa_R)),
        by = .(months_since_index, index_visit)][
          order(months_since_index)
        ][
          # Retain only first 7 years of followup
          months_since_index <=84 ] 
# fwrite(patients_7yr_exact, find_rstudio_root_file("data", "output", "7yr-followup-exact.csv"))




