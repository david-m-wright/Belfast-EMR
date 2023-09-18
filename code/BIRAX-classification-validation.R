# BIRAX AMD classification algorithm validation

library(tidyverse)

source(rprojroot::find_rstudio_root_file("code/BIRAX-assemble-cohort.R"))


# Identify two candidate groups of scans

# These selections are separate from fluid_history because they should start prior to the first injection 
# fluid_history is concerned only with treatment response

## Wet: "The first scan before the injection and also fluid appears in Notal's TRF column (both should apply at the same time)"
## Dry: "Scans before the scan that precedes the injection without fluid in Notal's TRF column."

# Details of all the examinations that were extracted
# Note that there were multiple ExamIDs and SeriesIDs (scans within exams) on some dates
oct_details <- fread(file = file.path(file_path, "OCT_ImageVariables.txt"))

# List of all dates on which an OCT volume scan was taken 
# When there are multiple exams on the same day, select the last (matches criteria used for raw NOA data)
# Assumes there was a problem with the earlier scan
oct_volumes <- unique(oct_details[ModalityType == "Volume" & ModalityProcedure == "IR_OCT" & ImageType == "OCT",], by = c("PatientID", "EyeCode", "ExamDate"), fromLast = TRUE)


## Generate a history of OCT scans for each eye 

# Add a variable to each data.table to make joins easier to navigate (otherwise one of the join columns will disappear)
# This is especially important for rolling joins
# Date format for the join columns must be identical
oct_volumes[, join_date := as.Date(ExamDate)]
eye_raw[,join_date := as.Date(index_date)]

# Set the join keys and sort the data - this must be done before determining OCT sequences
setkey(oct_volumes, PatientID, EyeCode, join_date)
setkey(eye_raw, PatientID, EyeCode, join_date)

# Mark each OCT scans with a sequence number, with zero being the baseline scan, taken on the day of injection
# Negative sequence numbers indicate scans taking place before baseline
oct_history <- eye_raw[oct_volumes, on = c("PatientID", "EyeCode")][,
  c("index_scan", 
    "oct_number") := list(ExamDate == index_date, 
                          seq_len(.N)), by = .(PatientID, EyeCode)][,
    oct_sequence := oct_number - oct_number[index_scan], by = .(PatientID, EyeCode)
  ]
oct_history[, history_date := as.Date(ExamDate)]
setkey(oct_history, PatientID, EyeCode, history_date)

# Example history
oct_history[
  ,  .(PatientID, EyeCode, join_date, index_date, ExamDate, history_date, index_scan, oct_sequence)
][
   PatientID == "0EC79835-C58D-CB89-8DF3-9C8B2C1C7AC4", # good example
  # PatientID == "0815B717-1564-D8C9-F606-3078FD3A74A3",
]


# Find baseline scans with evidence of fluid (>50nL) that have been successfully processed by the NOA
# Setup join
# fluid_raw[, join_date := as.Date(OCTDate)]
setkey(fluid_raw, PatientID, EyeCode, OCTDate)

# Make the join - matches on PatientID, EyeCode and date (not a rolling join) 
wet <- fluid_raw[oct_history, on = c("PatientID", "EyeCode", OCTDate ="history_date")][
  # Filter only scans at baseline and with fluid 
  oct_sequence == 0 & TRFVolumeNl > 50,][
    # A new column labelling these as wet
    , classif := "wet"
  ][,
    .(PatientID, EyeCode, ExamDate, ExamID, SeriesID, TRFVolumeNl, classif)
  ] 




## Finding dry scans - define these as scans that have taken place before baseline wet scans and have no evidence of fluid
# List of scans [xx] scans before the baseline wet scan to process in NOA (to find if dry)
avi_to_process <- oct_history[wet, on = c("PatientID", "EyeCode")][
  oct_sequence == -3, ][,
    .(PatientID, EyeCode, ExamDate, oct_sequence),][,
    # Generate file name stub to search for among unprocessed AVIs
  FileNameAVI := paste(PatientID, 
                           if_else(EyeCode == "R", "OD", "OS"),
                           format(ExamDate, format = "%Y%m%d"),
                           sep = "_")]

# List of AVI files that have not been processed
# For some patients there will be multiple scans on the same day
avi_available <- data.table(FileName=list.files(path = "F:\\Excluded_OCT_Scans", pattern = ".avi"))[, 
                  FileNameAVI := str_remove(FileName, "(_[0-9]{3})*.avi$")]
  
# Assemble priority list to process
avi_to_process <- avi_available[avi_to_process, on = "FileNameAVI"][
  !is.na(FileName), 
]

# # Move these to the new folder for NOA processing
 # fs::file_move(paste0("F:\\Excluded_OCT_Scans\\", avi_to_process$FileName), new_path = "F:\\RT-batch3\\HEYEX_Outputs")

# RT-batch1 is -1 scans (i.e. immediately before baseline)
# RT-batch2 is -2 scans



# Find the dry scans
# Make the join - matches on PatientID, EyeCode and date (not a rolling join) 
dry <- fluid_raw[oct_history, on = c("PatientID", "EyeCode", OCTDate ="history_date")][
  # Filter only scans at baseline and with fluid 
  oct_sequence == -2 & TRFVolumeNl == 0,][
    # A new column labelling these as wet
    , classif := "dry"
  ][,
    .(PatientID, EyeCode, ExamDate, ExamID, SeriesID, TRFVolumeNl, classif)
  ] 


# Generate list of the paired dry and wet scans
 dry_to_wet <-  data.table(bind_rows(
wet[dry, on = c("PatientID", "EyeCode"), nomatch = 0][,
  .(PatientID, EyeCode, ExamID, SeriesID, classif)
] ,
dry[wet, on = c("PatientID", "EyeCode"), nomatch = 0][,
                                         .(PatientID, EyeCode, ExamID, SeriesID, classif)
] ) )
 
 
 # Only include one eye from each patient
 dry_to_wet <- dry_to_wet[dry_to_wet[, row_number := seq_len(.N), by = .(PatientID)]
                          [row_number == 1, .(PatientID, EyeCode)], on = c("PatientID", "EyeCode")]
 

## Extract the .bmp files for scans

# Structure will be: classification (wet/dry) / PatientID 
# Naming convention? - keep as is, exclude the first frame of each OCT as this is the enface image

# Search needs to be on PatientID, EyeCode and ExamID and SeriesID


# Function to find all of the bmp files in the 'Volume' folder for a given PatientID
# Args: patient_id = patient ID whose folder to search
find_candidate_bmps <- function(patient_id){
  file.path("F:", "BIRAX_ProcessedOutputs", patient_id, "Volume") %>% 
    list.files() %>% 
    enframe(name = NULL, value = "file_name") %>% 
    mutate(PatientID = patient_id) %>% 
    separate(file_name, into = c("ExamID", "SeriesID", "EyeCode", "ModalityProcedure", "OCT", "ExamDate", "Frame", "Hash", "Extension"), remove = FALSE) %>% 
    mutate(across(Frame, as.numeric)) %>% 
    group_by(PatientID, EyeCode, ExamID, SeriesID) %>% 
    mutate(max_frames = max(Frame)-1) %>% 
    ungroup() %>% 
    as.data.table()
}


# Find list of candidate bmp file for each of the patients 
# Call unique() because both eyes may be needed form the same patient so avoid duplicating search
extracted_bmps <- vector("list", length = length(unique(dry_to_wet$PatientID)))

pb <- progress::progress_bar$new(
  format = "  Finding patients [:bar] :percent eta: :eta",
  total = length(extracted_bmps), clear = FALSE, width= 60)

for(i in 1:length(extracted_bmps)){
  pb$tick()
  extracted_bmps[[i]] <- find_candidate_bmps(unique(dry_to_wet$PatientID)[i]) 
}

candidate_bmps <- bind_rows(extracted_bmps)
candidate_bmps[, c("ExamID", "SeriesID") := lapply(.SD, as.numeric), .SDcols = c("ExamID", "SeriesID")]

# Match the scans by patient id, eye, ExamID and SeriesID
setkey(candidate_bmps, PatientID, EyeCode, ExamID, SeriesID)
setkey(dry_to_wet, PatientID, EyeCode, ExamID, SeriesID)

# If there are multiple scans at a single exam (i.e. identical PatientID, eye, date, number of frames),
# Selects the first in the series (lowest SeriesID)
# This does not match the selection criteria applied to the raw NOA data (would go for the most recent scan)
dry_to_wet_matched_bmps <- candidate_bmps[dry_to_wet][
  # , .SD[SeriesID == max(SeriesID)], by = ExamID][
    # Exclude the enface images (Frame 1 for each scan)
    Frame != 1
  ][,
    eye_id := paste(PatientID, if_else(EyeCode == "L", "OS", "OD"), sep = "_")
  ]


# Setup the directory structure for the bmp files
fs::dir_create(file.path("D:", "BIRAX-classification", "wet", unique(dry_to_wet_matched_bmps$eye_id)))
fs::dir_create(file.path("D:", "BIRAX-classification", "dry", unique(dry_to_wet_matched_bmps$eye_id)))

# Copy the bmp files into the new directories
fs::file_copy(file.path("F:", "BIRAX_ProcessedOutputs", dry_to_wet_matched_bmps$PatientID, "Volume", dry_to_wet_matched_bmps$file_name), 
              new_path = file.path("D:", "BIRAX-classification", dry_to_wet_matched_bmps$classif, dry_to_wet_matched_bmps$eye_id))




