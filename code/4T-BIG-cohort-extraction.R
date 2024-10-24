# Extraction of 4T-BIG cohort

source(rprojroot::find_rstudio_root_file("code/4T-BIG-cohort-construction.R"))



### Extract full imaging dataset for the selected cohort ###

# All injections as data.table for meta data
injections_dt <- data.table(injections_clean)[, .(PatientID, EyeCode, ExamDate = EncounterDate, injection = 1)]

# OCT scan metadata
oct_volume_metadata <- injections_dt[
  oct_volumes[patients_3yr_exact_list[, .(PatientID, EncounterDate, visit_sequence, months_since_index)], on = c("PatientID","ExamDate" = "EncounterDate"), nomatch = 0],
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


# out_dir <- find_rstudio_root_file("data", "output", "4T-BIG")
out_dir <- "E:\\4T-BIG"
# 
# 
# ## Copy imaging to the target directory, maintaining the source directory structure
# 
# # Create a directory for each individual
dir_create(paste0(out_dir, "\\Imaging\\", unique(patients_3yr_exact_list$PatientID)))
# 
# # Within each patient, create a directory for the OCT slices
dir_create(paste0(out_dir, "\\Imaging", unique(fnames$patientpath), "\\Volume"))
# 
# # Where present, create a directory for the multicolour enface images
oct_multicol_metadata <- oct_multicol[patients_3yr_exact_list[, .(PatientID, EncounterDate, visit_sequence, months_since_index)], on = c("PatientID","ExamDate" = "EncounterDate"), nomatch = 0]

oct_multicol_metadata[]
as.data.table(injections_clean)[, .(PatientID, EyeCode, EncounterDate)]

# [, injected := 1]

dir_create(paste0(out_dir, "\\Imaging", unique(str_extract(oct_multicol_metadata$FilePath,  ".*\\\\Single"))))
# 
# # Where present (and not already created for multicolour), create a directory for FA enface images 
# # (All are numbered slice 1)
oct_fa_metadata <- oct_fa[patients_3yr_exact_list[, .(PatientID, EncounterDate, visit_sequence, months_since_index)], on = c("PatientID","ExamDate" = "EncounterDate"), nomatch = 0]



dir_create(paste0(out_dir, "\\Imaging", unique(str_extract(oct_fa_metadata$FilePath,  ".*\\\\Single"))))
# 
# 
# # Populate with filtered scans
# # Multicolour (8 mins - 12Gb)
file_copy(paste0("F:\\BIRAX_ProcessedOutputs", oct_multicol_metadata$FilePath), paste0(out_dir, "\\Imaging", oct_multicol_metadata$FilePath))
# # FA
file_copy(paste0("F:\\BIRAX_ProcessedOutputs", oct_fa_metadata$FilePath), paste0(out_dir, "\\Imaging", oct_fa_metadata$FilePath))
# # OCT
file_copy(paste0("F:\\BIRAX_ProcessedOutputs", oct_slices$FilePath), paste0(out_dir, "\\Imaging", oct_slices$FilePath))
# 
# 
# 
# ## Output the metadata 
# 
dir_create(paste0(out_dir, "\\Metadata"))
# # For OCT, output a row for each slice.
oct_details[oct_slices[oct_volume_metadata[, .(volume_id, visit_sequence, months_since_index, injection)], on = "volume_id"], on = "FilePath", nomatch = 0] %>% 
fwrite(paste0(out_dir, "\\Metadata", "\\OCT-slice-metadata.csv"))
# # # For multiocolour and FA, a single row for each enface image
injections_dt[oct_multicol_metadata, on = c("PatientID", "EyeCode", "ExamDate")][, injection := if_else(is.na(injection), 0, injection)] %>% 
fwrite(paste0(out_dir, "\\Metadata", "\\multicol-metadata.csv"))
injections_dt[oct_fa_metadata, on = c("PatientID", "EyeCode", "ExamDate")][, injection := if_else(is.na(injection), 0, injection)] %>% 
fwrite(paste0(out_dir, "\\Metadata", "\\fa-metadata.csv"))
# # 
# # 14:53 multicol
# # 15:11 OCT
# 

# Eyes at baseline
oct_volume_metadata %>% filter(visit_sequence == 0, injection == 1)
oct_details[oct_slices[oct_volume_metadata[, .(volume_id, visit_sequence, months_since_index, injection)], on = "volume_id"], on = "FilePath", nomatch = 0] %>% 
  filter(visit_sequence == 0, injection == 1) %>% 
  count(PatientID, EyeCode)


# Patients at baseline
oct_volume_metadata %>% filter(visit_sequence == 0, injection == 1) %>% count(PatientID)
oct_details[oct_slices[oct_volume_metadata[, .(volume_id, visit_sequence, months_since_index, injection)], on = "volume_id"], on = "FilePath", nomatch = 0] %>% 
  filter(visit_sequence == 0, injection == 1) %>% 
  count(PatientID)



oct_volume_metadata %>% filter(visit_sequence == 0) %>% count(PatientID) %>% count(n)

oct_volume_metadata %>% filter(visit_sequence == 0) %>% count(PatientID, EyeCode, injection) %>% count(injection)


### Extract clinical data for the selected cohort ###



