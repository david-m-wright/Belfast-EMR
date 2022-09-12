# Assemble NOA dataset

library(tidyverse)
library(lubridate)
library(data.table)
library(rprojroot)

file_path <- find_rstudio_root_file("data/BIRAX data") 

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



### Sequencing of AVI files ###
# 
# # Generate list of AVI files to prioritise
# oct_to_process <- oct_visits %>% 
#   filter(is.na(`Analysis eligibility`)) %>% 
#   select(PatientID, EyeCode, EncounterDate) %>% 
#   # Generate the root filenames to search for
#   mutate(GenFileName = paste(PatientID, 
#                              if_else(EyeCode == "R", "OD", "OS"), 
#                              format(EncounterDate, format = "%Y%m%d"),       
#                              sep = "_"))


# write_csv(oct_to_process %>% 
#             select(GenFileName),
#           file = find_rstudio_root_file("data", "BIRAX data", "NOA-output", "oct_to_process.csv"))

### Prioritise OCTs to process ###
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

### Transfer half of OCTs for processing on second machine

# # All OCTs in the current batch (both processed and awaiting processing)
# oct_available <- list.files(path = "F:\\BIRAX\\HEYEX_Outputs", pattern = ".avi") %>%
#   enframe(name = NULL, value = "FileName") %>%
#   mutate(GenFileName = str_remove(FileName, ".avi$"))
# 
# # All OCTs processed (both succesfull and unsuccessful)
# oct_processed <- read_csv("F:\\BIRAX\\Final_19-07-2022-11-13-06\\Results.csv") %>% 
#   select(FileName) %>% 
#   full_join(read_csv("F:\\BIRAX\\Final_19-07-2022-11-13-06\\Logs\\Unsuccessful_scans.csv") %>% 
#               select(FileName),
#             by = "FileName")
# 
# # OCTs remaining to be processed
# oct_to_process <- anti_join(oct_available, oct_processed,
#           by = c("GenFileName" = "FileName")) %>% 
#    mutate(FileName = paste0("F:\\BIRAX\\HEYEX_Outputs\\", FileName)) 
# 
# # OCTs to move (half of those awaiting processing)
# oct_to_move <- oct_to_process %>% 
#   slice_tail(prop = 0.5)
# 
#  fs::file_move(pull(oct_to_move, FileName), new_path = "E:\\BIRAX\\HEYEX_Outputs")

### End of NOA machine section ###
