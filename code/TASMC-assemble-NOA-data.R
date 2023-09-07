# Assemble NOA dataset - TASMC
library(tidyverse)
library(lubridate)
library(data.table)
library(rprojroot)

# NOA data dictionary
noa_dictionary <- fread(find_rstudio_root_file("data-dictionary/NOA_dictionary.csv")) %>% 
  mutate(abbr_metric = str_squish(str_remove_all(Parameter, "Fluid|CST|Central|Sub Field|(T|t)emporal|Tempo|TS|Nasal|NS|(S|s)uperior|SS|(I|i)nferior|IS|Inner|Outer|3|6|\\(|\\)"))) 

# Assign parameter groups to character vectors for later use
noa_general <-  noa_dictionary %>% 
  filter(Category == "General", 
         !`Raw Parameter` %in% c("FileName", "Analysis eligibility")) %>% 
  pull(`Raw Parameter`)

noa_fluid <- noa_dictionary %>% 
  filter(Category == "Fluid Parameters",
         str_detect(`Raw Parameter`, "highest", negate = TRUE)) %>% 
  pull(`Raw Parameter`)

noa_retinal <- noa_dictionary %>% 
  filter(Category == "Retinal parameters"
         , str_detect(`Raw Parameter`, "highest|evidence", negate = TRUE)
  ) %>% 
  pull(`Raw Parameter`)

noa_grid <- noa_dictionary %>% 
  filter(Category == "ETDRS Grid Parameters") %>% 
  pull(`Raw Parameter`)


# Output from the Notal Ophthalmic Analyzer (NOA)

noa_raw <- bind_rows(
  # Patients with possible comordbidities visible on the OCTs (e.g. vein occlusion)
  # read_csv("//fas8200main-n2/OphBelfast/Results Notal_with original IDs_and_anon_185.csv") %>% 
  read_csv("//fas8200main-n2/OphBelfast/Final for analysis/Notal+MDC_2023-08-31_filtered_id_anon.csv") %>%   
    # select(-ID) %>% 
    # mutate(comorbidity = TRUE),
  # Patients with no comorbidities
  # read_csv("//fas8200main-n2/OphBelfast/Notal+MDC with additional IDs 16062023_anon_david.csv") %>% 
  # mutate(comorbidity = FALSE) %>% 
  rename_with(.fn = ~str_replace(., "\\(um\\)", "[um]")) %>%
  rename_with(.fn = ~str_replace(., "\\(nl\\)", "[nl]")) %>% 
  rename_with(.fn = ~str_replace(., "\\(mm\\^2\\)", "[mm\\^2]"))
  ) %>% 
  rename(PatientID = ID_anon, 
         EyeCode = Laterality) %>% 
  mutate(across(DATE, ymd),
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
# Selected scans were eligible
# Fewest missing values (indicator of neural network processing),
# Highest number of frames in the AVI (B-scans)
# If tied, the final scan taken on that date is selected. 
noa <- noa_raw %>% 
  filter(`Analysis eligibility` == 1) %>% 
  filter(FramesNum >= 18) %>% 
  mutate(OCTMissing = rowSums(across(!matches("Bscan"), is.na))) %>% 
  group_by(PatientID, EyeCode, DATE) %>% 
  slice_max(`Analysis eligibility`, n = 1, with_ties = TRUE) %>% 
  slice_min(OCTMissing, n = 1, with_ties = TRUE) %>% 
  slice_max(FramesNum, n = 1, with_ties = TRUE) %>% 
  slice_tail(n = 1) %>% 
  ungroup()

