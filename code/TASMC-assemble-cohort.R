# Setup TASMC cohort
library(tidyverse)
library(data.table)
library(conflicted)

conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")


source(rprojroot::find_rstudio_root_file("code/EMR-helper-functions.R"))
source(rprojroot::find_rstudio_root_file("code/TASMC-assemble-NOA-data.R"))

# Derived variable labels
var_desc <- read_csv(find_rstudio_root_file("data-dictionary/Belfast-EMR-derived-variables.csv"))

# Load the clinical data (MD Clone)
mdc_raw <- read_csv("//fas8200main-n2/OphBelfast/MDC_Updated_16062023_anon_david.csv") %>%
  rename(PatientID = ID_anon, 
         EyeCode = EYE,
         Gender = GENDER) %>% 
  mutate(Date = dmy(Date),
         exclude_research_patient = is.na(`Is_RP_Flg RP= Research patients 1=Yes 0=No`)|  
           `Is_RP_Flg RP= Research patients 1=Yes 0=No` == 1,
         exclude_covid = `Is_Pre_COVID_Flg  1=Yes 0=No` == 0) %>% 
  mutate(across(Gender, as.factor)) %>% 
  filter(!is.na(Date),
         !is.na(AGE))

# Exclude research patients and those first treated during Covid
mdc <- mdc_raw %>% 
  filter(!exclude_research_patient,
         !exclude_covid)

# Patient details
patients_raw <- mdc %>% 
  distinct(PatientID, Gender)


# Find date of first injection for each eye (index_date) and calculate treatment intervals
injections_raw <- mdc %>% 
  filter(EVENT %in% c("AVASTIN", "EYLEA", "LUCENTIS")) %>% 
  group_by(PatientID, EyeCode) %>% 
  mutate(index_date = min(Date),
         years_treated = interval(index_date, Date)/dyears(),
         # Year of treatment for each eye starts at 1
         treatment_year = as.factor(floor(years_treated) + 1)) %>% 
  ungroup() %>% 
  select(PatientID, EyeCode, AGE, index_date, Date, EVENT, years_treated, treatment_year)


## Prepare series of measurements over time ##

# Injection history summary by eye
# Find date of final injections and number of injections for each eye
injection_summary_eye <- injections_raw %>% 
  group_by(PatientID, EyeCode) %>% 
  summarise(index_date = min(Date),
            final_injection_date = max(Date),
            total_injections = n(), 
            .groups = "drop") %>% 
  mutate(total_intervals = total_injections - 1)


# Visual acuity
visual_acuity <- mdc %>%
  filter(EVENT == "VA") %>%
  select(PatientID, Date, `BCVA OD`, `BCVA OS`) %>%
  pivot_longer(
    cols = c("BCVA OD", "BCVA OS"),
    names_prefix = "BCVA ",
    names_to = "EyeCode",
    values_to = "va_logmar"
  ) %>%
  group_by(PatientID, EyeCode, Date) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(
    va_category_snellen = cut(
        va_logmar,
        breaks = c(min(va_logmar, na.rm = TRUE),
          0.29,
          0.59,
          0.99,
          max(va_logmar, na.rm = TRUE)
        ),
        labels = c(
          "Good (VA>6/12)",
          "Moderate (6/24 < VA <= 6/12)",
          "Partially sighted (6/60 < VA <= 6/24)",
          "Blind (VA <= 6/60)"
        ),
        include.lowest = TRUE)
  )
 

# For each visual acuity measurement, calculate the time since the first injection (index_date)
va_raw <- as.data.table(injection_summary_eye)[, .(PatientID, EyeCode, index_date)][ 
  # Join injection_summary_eye and visual_acuity tables
  as.data.table(visual_acuity), .(PatientID, EyeCode, index_date, Date, va_logmar, va_category_snellen), on = .(PatientID, EyeCode)][
    # Calculate months since index date
    , months_since_index:=interval(index_date, Date)/dmonths()][
      # Calculate years ince index date
      , years_since_index := interval(index_date, Date)/dyears()][  
        # Mark baseline measurements (closest measurement to baseline)  
         , baseline := abs(months_since_index) == min(abs(months_since_index)), by = .(PatientID, EyeCode)][
           # Only retain eyes that received injections (and hence had a baseline measurement)
           !is.na(baseline),
           ]


## NOA fluid measurements
# For each OCT fluid measurement, calculate the time since the first injection (index_date)
fluid_raw <- as.data.table(injection_summary_eye)[
  , .(PatientID, EyeCode, index_date)][
    # Join injection_summary_eye and NOA fluid measurements
    as.data.table(noa), on = .(PatientID, EyeCode)][
      # Calculate months since index date
      , months_since_index:=interval(index_date, DATE)/dmonths()][
        # Mark baseline measurements (closest measurement to baseline)  
        , baseline := abs(months_since_index) == min(abs(months_since_index)), by = .(PatientID, EyeCode)]


### Assemble eye level dataset (single row per eye) ###

eye_raw <- patients_raw %>% 
  
  # Extract age at first injection (index date)
  inner_join(injections_raw %>% 
              filter(Date == index_date) %>% 
            select(PatientID, EyeCode, index_age = AGE),
            by = "PatientID", multiple = "all") %>% 
  
  # Calculate years treated (index to final injection)
  left_join(injection_summary_eye, by = c("PatientID", "EyeCode")) %>% 
  mutate(years_treated = interval(index_date, final_injection_date)/dyears()) %>% 
  
  # Calculate years from index to end of VA monitoring
  left_join(va_raw[!is.na(years_since_index), by = .(PatientID, EyeCode), .(years_va = max(years_since_index))],
            by = c("PatientID", "EyeCode")) %>% 
  
  # For a proportion of eyes years_treated is greater than years_observed 
  # Where this occurs, set years observed to the greater of the two
  mutate(years_observed = pmax(years_treated, years_va)) %>% 
  
  # Join VA measured on the index date or up to 40 days prior as the baseline
  left_join(va_raw %>%
              # filter(baseline) %>%
               filter(baseline, months_since_index > -1.3 & months_since_index <=0) %>%
              select(PatientID, EyeCode, va_logmar, va_category_snellen),
            by = c("PatientID", "EyeCode")) %>% 
  
  
  # Join presence of fluid measurements (closest to index date but up to 40 days prior)
  left_join(fluid_raw %>% 
               filter(baseline) %>% 
               # filter(baseline, months_since_index > -1.3 & months_since_index <=0) %>% 
              transmute(PatientID, EyeCode, fluid_measurement = baseline), 
            by = c("PatientID", "EyeCode")) %>% 
  
  mutate(across(fluid_measurement, ~if_else(is.na(.), FALSE, .))) %>%
  
  # Exclusion criteria
  mutate(
    
    # 50 or older on date of first injection (index date)
    exclude_age = index_age < 50,
    
    # Exclude <3 anti-VEGF injections
    exclude_lt3_injections = if_else(total_injections < 3 | is.na(total_injections), TRUE, FALSE),
    
    # No injections after index date (so no treatment intervals can be calculated)
    exclude_no_intervals = total_intervals == 0,
    
    # Those without a baseline VA measurement
    exclude_no_va = is.na(va_logmar),
    
    # Those without a baseline fluid measurement
    exclude_no_fluid = !fluid_measurement
  ) %>% 
  
  # Add an eye id to allow joining of model predictions with eye attibutes
  mutate(eye_id = factor(paste(PatientID, EyeCode, sep = "-")))


# Check exclusions
eye_raw %>% 
  count(across(matches("exclude")))
eye_raw %>% count(exclude_no_va)


# Apply the exclusion criteria
# Retain those with no baseline VA for use in other analyses
eye <- eye_raw %>% 
  filter(!exclude_age, 
         !exclude_lt3_injections, 
         !exclude_no_intervals, 
         # !exclude_no_va, 
         !exclude_no_fluid) %>% 
  select(-matches("exclude")) %>% 
  
  # Derived variable categorising index age into decades
  mutate(Age = cut(index_age, right = F, include.lowest = T, breaks = c(50, 60, 70, 80, 90, max(index_age))),
         Age = IntervalToInequality(Age, "yrs"))



## Prepare time series of injections, VA, OCT and fluid measurements for the selected eyes ##

# Just the injection history for the selected eyes
injections <- injections_raw %>% 
  inner_join(eye %>% 
               select(PatientID, EyeCode), by = c("PatientID", "EyeCode"))

# Just the VA history relating to the selected eyes
va_history_raw <- va_raw %>% 
  inner_join(eye %>% 
               select(PatientID, EyeCode), by = c("PatientID", "EyeCode")) %>% 
  # Roll forward baseline VA measurements so that all series start at the index date 
  # (max roll forward is 40 days)
  # Would gain ~20 eyes if allowed two week roll back after index date.
  mutate(months_since_index = if_else(baseline & months_since_index > -1.3 & months_since_index < 0, 0, months_since_index)) %>% 
  # Exclude all measurements prior to baseline
  filter(months_since_index >= 0) %>% 
  mutate(treatment_months = ceiling(months_since_index))

# Take snapshots of VA at various stages of follow-up
# Find closest VA measurement to each of the specified snapshot times, within set tolerances
# Follow-up times in months at which the VA status is to be taken, along with acceptable tolerances
# Baseline dates must match exactly
snapshots <- data.table(months_since_index = c(0, 6, 12, 24, 36, 48, 60, 84, 120), 
                        follow_up_month = c(0, 6, 12, 24, 36, 48, 60, 84, 120), 
                        tolerance = c(0, 2, 2, 2, 2, 2, 2,2,2),
                        key = "months_since_index") %>% 
  .[, `:=` (lower_lim = follow_up_month - tolerance,
            upper_lim = follow_up_month + tolerance,
            yr = follow_up_month/12,
            follow_up_year = as.factor(follow_up_month/12))]


# Set keys to join on (so the .on() specification is not needed in the join)
setkey(snapshots, months_since_index)
setkey(va_history_raw, months_since_index)

# Identify the VA measurement closest to each of a given set of snapshots for each eye
# Retain all measurements but mark those that can be used for a snapshot
# Candidate VA measurement must be within a given tolerance either side of the snapshot time
# Break ties of equal distance from the snapshot by choosing the earlier VA measurement (rolling back)
va_history <- snapshots[va_history_raw, roll = "nearest"][
  , snapshot := .I == .I[which.min(abs(months_since_index - follow_up_month))] &
    months_since_index >= lower_lim &
    months_since_index <= upper_lim, by = c("PatientID", "EyeCode", "follow_up_month")][
      
      # Calculate change in VA from baseline
      , `:=` (va_change_logmar = va_logmar - va_logmar[(snapshot & baseline)])
      , by=.(PatientID, EyeCode)][
        # Classify VA change into three categories
        , va_change_lines := cut(va_change_logmar, 
                                 breaks = c(min(va_change_logmar, na.rm = TRUE), -0.2, 0.19, max(va_change_logmar, na.rm = TRUE)), 
                                 labels = c("gained >=2 lines", "< 2 lines change", "lost >= 2 lines"), 
                                 include.lowest = TRUE)][
                                   # Add the treatment year to match up with injection counts
                                   , treatment_year := as.factor(floor(years_since_index)+1)]



# Fluid history for each eye
fluid_history_raw <- fluid_raw %>% 
  inner_join(eye %>% 
               filter(fluid_measurement) %>% 
               select(PatientID, EyeCode), by = c("PatientID", "EyeCode")) %>% 
  # Roll forward baseline fluid measurements up to 40 days prior to baseline
  mutate(months_since_index = if_else(baseline & months_since_index > -1.3 & months_since_index < 0, 0, months_since_index)) %>%
  # Assign each measurement to a treatment month (baseline = zero, within first month of treatment = 1)
  mutate(treatment_months = ceiling(months_since_index)) %>% 
  # Exclude all measurements prior to baseline
   filter(months_since_index >= 0) %>% 
  
  # Classification of AMD type by presence of SRF and IRF
  mutate(amd_type = factor(case_when(SRF == 1 & IRF == 1 ~ "SRF and IRF",
                                     SRF == 0 & IRF == 1 ~ "IRF only",
                                     SRF == 1 & IRF == 0 ~ "SRF only",
                                     SRF == 0 & IRF == 0 ~ "No fluid"), levels = c("No fluid", "SRF only", "IRF only", "SRF and IRF")))


# Identify thickness measurements corresponding to snapshots
setkey(fluid_history_raw, months_since_index)
fluid_history <- snapshots[fluid_history_raw, roll = "nearest"][
  , snapshot := .I == .I[which.min(abs(months_since_index - follow_up_month))] &
    months_since_index >= lower_lim &
    months_since_index <= upper_lim, by = c("PatientID", "EyeCode", "follow_up_month")]#[


# # Output the fluid history for exploratory analysis
# fluid_output <- fluid_history %>%
#   filter(baseline) %>%
#   inner_join(eye, by = c("PatientID", "EyeCode", "index_date")) %>% 
#   rename(`First_Inj_Date` = `1st_Inj_Date`) 
# names(fluid_output) <- str_replace_all(names(fluid_output), "(\\[|\\]|\\(|\\)|[[:space:]])|\\?", "_") 
# names(fluid_output) <- str_replace_all(names(fluid_output), "\\^|\\=", "")   
# # write_csv(fluid_output, "//fas8200main-n2/OphBelfast/fluid-history.csv")
# haven::write_sav(fluid_output, "//fas8200main-n2/OphBelfast/fluid-history.sav")



treatment_intervals <- injections %>%
  group_by(PatientID, EyeCode) %>% 
  mutate(next_encounter = lead(Date),
         treatment_interval_days = interval(Date, next_encounter)/ddays(), 
         treatment_interval_weeks = interval(Date, next_encounter)/dweeks(),
         # Indicator of short treatment interval 
         treatment_interval_short = treatment_interval_days <= 40) %>% 
  filter(!is.na(next_encounter)) %>% 
  # Timeline since the index date for each person
  # First and second year since treatment start concept not relevant because do not know the date of first injection (no clinic entry date recorded)
  mutate("Treatment interval" = if_else(treatment_interval_weeks > 12, ">12 weeks", "<=12 weeks"),
         # Calculate month and year of treatment for each eye (since index date)       
         treatment_months = interval(index_date, Date)/dmonths(),
         treatment_year = as.factor(floor(interval(index_date, Date)/dyears()) + 1)) %>% 
  # Injection sequence (index injection = 1)
  mutate(injection_seq = row_number()) %>% 
  ungroup() %>% 
  mutate(sequence_id = abbreviate(paste(PatientID, EyeCode))) 


## Compare visual acuity and fluid measurements at baseline
va_fl <- eye %>% 
  inner_join(fluid_history %>% 
               filter(baseline),
             by = c("PatientID", "EyeCode"))

va_fl_correlations <- va_fl %>% 
  select(all_of(c(noa_fluid, noa_retinal, noa_grid)) & where(is.numeric)) %>% 
  map(~cor.test(.x, va_fl$va_logmar)) %>% 
  map_dfr(~broom::tidy(.), .id = "Variable") %>% 
  mutate(across(p.value, ~format.pval(., eps = 0.001, digits = 2)),
         across(c(estimate, conf.low, conf.high), .fns = ~formatC(.x, format = "f", digits = 2))) %>% 
  mutate(`Correlation with VA` = paste0(estimate, " (", conf.low, ", ", conf.high, ")")) %>% 
  transmute(Variable, estimate, `Correlation with VA`,  P = p.value)

# Correlations between VA and fluid measurements only among those with IRF AND SRF at baseline
va_fl_correlations_srf_irf <- va_fl %>% 
  filter(amd_type == "SRF and IRF") %>% 
  select(all_of(c(noa_fluid, noa_retinal, noa_grid)) & where(is.numeric)) %>% 
  map(~cor.test(.x, filter(va_fl, amd_type == "SRF and IRF")$va_logmar)) %>% 
  map_dfr(~broom::tidy(.), .id = "Variable") %>% 
  mutate(across(p.value, ~format.pval(., eps = 0.001, digits = 2)),
         across(c(estimate, conf.low, conf.high), .fns = ~formatC(.x, format = "f", digits = 2))) %>% 
  mutate(`Correlation with VA` = paste0(estimate, " (", conf.low, ", ", conf.high, ")")) %>% 
  transmute(Variable, estimate, `Correlation with VA`,  P = p.value)


# Fluid measurements at baseline by ETDRS grid in long format
baseline_fluid_grid <- fluid_history %>% 
  filter(baseline) %>% 
  select(all_of(noa_grid)) %>% 
  pivot_longer(cols = everything()) %>% 
  group_by(name) %>% 
  summarise(valform = paste0(round(mean(value, na.rm = TRUE), 1), " \n(", round(sd(value, na.rm = TRUE), 1), ")")) %>% 
  mutate(region = toupper(str_sub(str_extract(name, "CS|(Tempo|Nasal|Superior|Inferior)(?=3|6)"), 1, 3)),
         circ = str_extract(name, "3|6"),
         circ = as.numeric(if_else(is.na(circ), "1", circ))) %>% 
  inner_join(noa_dictionary %>% 
               filter(Category == "ETDRS Grid Parameters"),
             by = c("name" = "Raw Parameter"))



