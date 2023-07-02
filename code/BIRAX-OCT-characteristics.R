# OCT scan characteristics from the EMR cohort

library(tidyverse)
library(lubridate)
# data.table for speed in lead/lag operations
library(data.table)
library(rprojroot)
library(janitor)
library(conflicted)
conflict_prefer("filter", "dplyr")


source(find_rstudio_root_file("code/EMR-helper-functions.R"))

# Load the OCT scan list (slow - 16M records)
oct_images <- fread("F:/BIRAX/OCT_ImageVariables.txt")

oct_images %>% 
  count(ModalityType, ModalityProcedure, ImageType, Resolution) %>% 
  arrange(desc(n)) %>% 
  adorn_totals() %>% 
  adorn_percentages("col") %>% 
  adorn_pct_formatting() %>% 
  adorn_ns("front")

# Any exclusions here?
oct_images %>% 
  + count(ExaminedStructure)

# Volumes or movies?
# volumes were more frequent


# Modalities taken by eye and visit date

# There are multiple files for a single scan (e.g. jpgs for each B-scan)
# Set keys to speed up sorts 
setkeyv(oct_images, c("PatientID", "EyeCode", "ExamDateTime", "ModalityProcedure", "Modality", "ModalityType"))

# Individual exams
unique(oct_images, by = c("PatientID", "EyeCode", "ExamDateTime")) 

# Modalities captured at each exam
dcast(oct_images, PatientID + EyeCode + ExamDateTime ~ModalityProcedure +  Modality + ModalityType, value.var = "ModalityProcedure", fun.aggregate = length)

# Cannot rely on ExamID
unique(oct_images, by = "ExamID")

# Just volumes and movies
vm <- oct_images[ModalityType %in% c("Volume", "Movie")] %>% 
  dcast(PatientID + EyeCode + ExamDateTime ~ModalityProcedure +  Modality + ModalityType, value.var = "ModalityProcedure", fun.aggregate = length) 

vm %>% 
  .[, c(-1, -2, -3)] %>% colSums()
  
vm[, lapply(.SD, function(x){sum(x>0)}), .SDcols = c("FA_OCT_OCT_Volume", "ICGA_OP_Movie")]

# Convert file counts to binary
vmb <- vm[, lapply(.SD, function(x){sum(x>0)}), by = .(PatientID, EyeCode, ExamDateTime)]

vmb %>% 
  .[, c(-1, -2, -3)] %>% 
  colSums()

vmb %>% 
  count(across(.cols = 4:ncol(.))) %>% 
  arrange(desc(n))


# Modalities taken
unique(oct_images, by = c("PatientID", "EyeCode", "ExamDateTime", "ModalityProcedure", "ModalityType")) %>% 
  count(ModalityProcedure, ModalityType)

oct_images[1, .(PatientID, EyeCode, ExamDateTime)][oct_images, on = .(PatientID, EyeCode, ExamDateTime), nomatch = NULL] %>% 
  count(Modality, ModalityProcedure, ModalityType, ImageType, LightSource, Resolution)

oct1 <- oct_images[1, .(PatientID, EyeCode, ExamDateTime)][oct_images, on = .(PatientID, EyeCode, ExamDateTime), nomatch = NULL] 

dcast(oct1, PatientID + EyeCode + ExamDateTime ~ModalityProcedure +  Modality + ModalityType, value.var = "ModalityProcedure", fun.aggregate = length)
  
count(oct1, Modality, ModalityProcedure, ModalityType)


# Multiple exams per day for any patient?
# Does not matter at this stage


oct_by_type <- oct_images[ModalityType == "Volume", .(PatientID, EyeCode, ExamDate, ModalityProcedure)]

dcast(oct_by_type, PatientID + EyeCode + ExamDate ~ ModalityProcedure, value.var = "ModalityProcedure")[1,][
  oct_images, on = .(PatientID, EyeCode, ExamDate), nomatch = NULL 
]
  # inner_join(oct_images, by = c("PatientID", "EyeCode", "ExamDate"))
# How to do the join with data.table?    



