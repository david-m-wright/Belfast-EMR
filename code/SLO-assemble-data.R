# Analysis of scanning laser ophthalmoscopy (SLO) images using SLOctolyzer (University of Edinburgh)

library(data.table)
library(stringr)
library(openxlsx)


# Prepare dataset

# Details of all the examinations that were extracted
# Note that there were multiple ExamIDs and SeriesIDs (scans within exams) on some dates
oct_details <- fread(file = "F:\\OCT_ImageVariables.txt")
oct_details[, c("StartCoordX", "StartCoordY") := lapply(.SD, as.numeric), .SDcols = c("StartCoordX", "StartCoordY")]
setkey(oct_details, FilePath)

# List of all dates on which an OCT volume scan was taken 
# When there are multiple exams on the same day, select the last (matches criteria used for raw NOA data)
# Assumes there was a problem with the earlier scan
oct_volumes <- unique(oct_details[ModalityType == "Volume" & ModalityProcedure == "IR_OCT" & ImageType == "OCT",], by = c("PatientID", "EyeCode", "ExamDate"), fromLast = TRUE)

# Find the en face images to go with each volume
oct_slo <- unique(oct_details[ModalityType == "Volume" & ModalityProcedure == "IR_OCT" & ImageID == 1,], by = c("PatientID", "EyeCode", "ExamDate"), fromLast = TRUE)
nrow(oct_volumes) == nrow(oct_slo)

# Copy all SLO files to a separate directory for processing
# fs::file_copy(paste0("F:\\BIRAX_ProcessedOutputs", oct_slo$FilePath), new_path = "D:/SLO")

# Generate metadata with lateralities and pixel scaling (do not have location)
oct_slo_meta <- 
  oct_slo[, c("Filename", "Scale", "Location", "Eye") := #list(str_remove(basename(FilePath), ".bmp"),
                                                              list(basename(FilePath),
                                                                   
                                                           ScaleX*1000,
                                                           "",
                                                           ifelse(EyeCode=="R", "Right", "Left"))][
                                                            , c("Filename", "Scale", "Location", "Eye")
                                                           ]
# write.xlsx(oct_slo_meta, "D:/SLO/fname_resolution_location_eye.xlsx", sheetName="fname_resolution_location_eye")

# Can run using reticulate or at the console (less likely to crash R)
# library(reticulate)
# use_condaenv("slo-analysis")
# reticulate::py_run_file("D:\\OneDrive - Queen's University Belfast\\Work\\SLOCT\\SLOctolyzer_demo_300524\\sloctolyzer\\main.py")
# Start run 04/06/2024 - estimated 30 day runtime



