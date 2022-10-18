# Assemble BIRAX analysis cohort (eye level)

library(tidyverse)
library(lubridate)

source(rprojroot::find_rstudio_root_file("code/BIRAX-assemble-Medisoft-data.R"))
source(rprojroot::find_rstudio_root_file("code/BIRAX-assemble-NOA-data.R"))


# Find first AMD diagnosis date for each eye
AMD_patients <- copath %>% 
  filter(CoPathologyCode %in% c("R316", "R337")) %>% 
  select(PatientID, EyeCode, DiagnosisDescription_AMD = CoPathologyDesc, DateofDiagnosis_AMD = EncounterDate) %>% 
  bind_rows(
    diagnoses %>% 
      filter(str_detect(DiagnosisDescription, "AMD|age-related macular degeneration"), 
             !str_detect(DiagnosisDescription, "suspected|no evidence|dry"),
             !is.na(DiagnosisICD10Code)) %>% 
      select(PatientID, EyeCode, DiagnosisDescription_AMD = DiagnosisDescription, DateofDiagnosis_AMD = DateofDiagnosis)) %>% 
  group_by(PatientID, EyeCode) %>% 
  slice_min(order_by =  DateofDiagnosis_AMD, with_ties = FALSE) %>% 
  ungroup()
# Note small number with dry AMD at first diagnosis.
# AMD_patients %>% count(DiagnosisDescription_AMD)

# Find RVO patients
RVO_patients <- copath %>% 
  filter(CoPathologyCode == "R318") %>% 
  select(PatientID, EyeCode, DiagnosisDescription_RVO = CoPathologyDesc, DateofDiagnosis_RVO = EncounterDate) %>% 
  bind_rows(
    diagnoses %>% 
      filter(str_detect(DiagnosisDescription, "occlusion")) %>% 
      select(PatientID, EyeCode, DiagnosisDescription_RVO = DiagnosisDescription, DateofDiagnosis_RVO = DateofDiagnosis)) %>% 
  group_by(PatientID, EyeCode) %>% 
  slice_min(order_by = DateofDiagnosis_RVO, with_ties = FALSE) %>% 
  ungroup()

# Find DMO/DR patients
DMO_DR_patients <- bind_rows(
  copath %>% 
    filter(CoPathologyCode %in% c("R318")) %>% 
    select(PatientID, EyeCode, DiagnosisDescription_DR = CoPathologyDesc, DateofDiagnosis_DR = EncounterDate),
  
  diagnoses %>% 
    filter(str_detect(DiagnosisDescription, "diabetic macular oedema"),
           !str_detect(DiagnosisDescription, "absent")) %>% 
    select(PatientID, EyeCode, DiagnosisDescription_DR = DiagnosisDescription, DateofDiagnosis_DR = DateofDiagnosis),
  
  dr_assessment %>% 
    select(PatientID, EyeCode, DiagnosisDescription_DR = DRAssessmentDesc, DateofDiagnosis_DR = EncounterDate) %>% 
    filter(DiagnosisDescription_DR != "none")
  
) %>% 
  group_by(PatientID, EyeCode) %>% 
  slice_min(order_by = DateofDiagnosis_DR, with_ties = FALSE) %>% 
  ungroup()



### Prepare series of measurements over time ###

# These analyses are all based on the injection history for each eye, starting at the date of the first injection

# Injection history summary by eye
# Find date of first injection (IndexDate) and final injection
# and number of injections for each eye
injection_summary_eye  <- injections_clean %>% 
  group_by(PatientID, EyeCode) %>% 
  summarise(index_date = min(EncounterDate),
            final_injection_date = max(EncounterDate), 
            total_injections = n(),
            .groups = "drop") %>% 
  mutate(total_intervals = total_injections - 1)


# For each visual acuity measurement, calculate the time since the first injection (index_date)
va_raw <- as.data.table(injection_summary_eye)[, .(PatientID, EyeCode, index_date)][ 
  # Join injection_summary_eye and visual_acuity tables
  as.data.table(visual_acuity), .(PatientID, EyeCode, index_date, EncounterDate, va_logmar, va_etdrs, va_category_snellen), on = .(PatientID, EyeCode)][
  # Calculate months since index date
  , months_since_index:=as.interval(index_date, EncounterDate)/dmonths()][
  # Mark baseline measurements (closest measurement to baseline)  
  , baseline := abs(months_since_index) == min(abs(months_since_index)), by = .(PatientID, EyeCode)]


## OCT thickness metrics
oct_thickness_raw <- as.data.table(injection_summary_eye)[
  , .(PatientID, EyeCode, index_date)][
  # Join injection_summary_eye and oct_thickness tables
  as.data.table(oct_thickness), on = .(PatientID, EyeCode)][
  # Calculate months since index date
  , months_since_index:=as.interval(index_date, ExamDate)/dmonths()][
  # Mark baseline measurements (closest measurement to baseline)  
  , baseline := abs(months_since_index) == min(abs(months_since_index)), by = .(PatientID, EyeCode)]

# Calculate changes in thickness measurements from baseline 
# Could do this with a set of baseline columns and a set of change columns
# ss <- oct_thickness_raw %>% 
#   .[, .(PatientID, EyeCode, ExamDate, OuterSuperior_Volume, OuterInferior_Volume)] 
#   cnames <- c("OuterSuperior_Volume", "OuterInferior_Volume")
# Change x[1] to x[baseline]
#   ss[, .SDcols = cnames,  (paste0(cnames, "_change")) := lapply(.SD, function(x) x - x[1]), by=list(PatientID, EyeCode)] %>% 
#     print()


## NOA fluid measurements
# For each OCT fluid measurement, calculate the time since the first injection (index_date)
fluid_raw <- as.data.table(injection_summary_eye)[
    , .(PatientID, EyeCode, index_date)][
    # Join injection_summary_eye and oct_thickness tables
    as.data.table(noa), on = .(PatientID, EyeCode)][
    # Calculate months since index date
    , months_since_index:=as.interval(index_date, OCTDate)/dmonths()][
    # Mark baseline measurements (closest measurement to baseline)  
    , baseline := abs(months_since_index) == min(abs(months_since_index)), by = .(PatientID, EyeCode)]


### Assemble eye level dataset (single row per eye) ###

eye_raw <- patients %>% 
  select(SiteID, PatientID, PerturbedDateofBirth, PerturbedCurrentAge, PerturbedAgeatDeath, Gender) %>% 
  
  # Add diagnoses for AMD, RVO and DMO/DR
  left_join(AMD_patients, by = "PatientID") %>% 
  left_join(RVO_patients, by = c("PatientID", "EyeCode")) %>% 
  left_join(DMO_DR_patients, by = c("PatientID", "EyeCode")) %>% 
  
  # Calculate age at first injection (index date)
  left_join(injection_summary_eye, by = c("PatientID", "EyeCode")) %>% 
  mutate(index_age =  as.interval(PerturbedDateofBirth, index_date)/dyears()) %>% 
  # Calculate years observed (index to current age/death)
  mutate(years_observed = coalesce(PerturbedCurrentAge, PerturbedAgeatDeath) - index_age) %>% 
  # Calculate years treated (index to final injection)
  mutate(years_treated = interval(index_date, final_injection_date)/dyears()) %>% 
  
  # Join VA measured on the index date or up to 2 weeks prior as the baseline
  left_join(va_raw %>%
              filter(baseline, months_since_index > -0.5 & months_since_index <=0) %>% 
              transmute(PatientID, EyeCode, va_logmar, "ETDRS letters" = va_etdrs),
            by = c("PatientID", "EyeCode")) %>% 
  
  # Join presence of OCT thickness measurements (closest to index date but up to 2 weeks prior)
  # Even though an eye may have no baseline (indicating not measured here) it may have later measurements that can be included
  # although change metrics cannot be calculated for these.
  left_join(oct_thickness_raw %>% 
              filter(baseline, months_since_index > -0.5 & months_since_index <=0) %>% 
              transmute(PatientID, EyeCode, thickness_measurement = baseline), 
            by = c("PatientID", "EyeCode")) %>% 
  
  # Join presence of fluid measurements (closest to index date but up to 2 weeks prior)
  left_join(fluid_raw %>% 
              filter(baseline, months_since_index > -0.5 & months_since_index <=0) %>% 
              transmute(PatientID, EyeCode, fluid_measurement = baseline), 
            by = c("PatientID", "EyeCode")) %>% 
  
  mutate(across(c(thickness_measurement, fluid_measurement), ~if_else(is.na(.), FALSE, .))) %>%
  
  # Exclusion criteria
  mutate(
    # Exclude no AMD diagnosis
    exclude_no_AMD = is.na(DiagnosisDescription_AMD),
    
    # 50 or older on date of first injection (index date) (or no injections)
    exclude_age = index_age < 50 | is.na(index_age),
    
    # Perturbation of ages makes follow-up period invalid (negative)
    exclude_age_perturb = years_observed < 0,
    
    # Exclude <3 anti-VEGF injections
    exclude_lt3_injections = if_else(total_injections < 3 | is.na(total_injections), TRUE, FALSE),
    
    # Diagnosis of DMO/DR or RVO
    exclude_DR_DMO = !is.na(DiagnosisDescription_DR),
    exclude_RVO = !is.na(DiagnosisDescription_RVO),
    
    # No injections after index date (so no treatment intervals can be calculated)
    exclude_no_intervals = total_intervals == 0,
    
    # Those without a baseline VA measurement
    exclude_no_va = is.na(va_logmar)) 


# Apply the exclusion criteria
eye <- eye_raw %>% 
  filter(!exclude_no_AMD, !exclude_age, !exclude_age_perturb, !exclude_lt3_injections,  !exclude_RVO, !exclude_DR_DMO, !exclude_no_va) %>% 
  select(-matches("exclude")) %>% 
  
  # Derived variable categorising index age into decades
  mutate(Age = cut(index_age, right = F, include.lowest = T, breaks = c(50, 60, 70, 80, 90, max(index_age))),
       Age = IntervalToInequality(Age, "yrs"))
  

## Prepare time series of VA, OCT and fluid measurements for the selected eyes ##

# Just the VA history relating to the selected eyes
va_history_raw <- va_raw %>% 
  inner_join(eye %>% 
               select(PatientID, EyeCode), by = c("PatientID", "EyeCode")) %>% 
  # Roll forward baseline VA measurements so that all series start at the index date 
  # (max roll forward is 2 weeks)
  mutate(months_since_index = if_else(baseline & months_since_index > -0.5 & months_since_index < 0, 0, months_since_index)) %>% 
  # Exclude all measurements prior to baseline
  filter(months_since_index >= 0)

# Take snapshots of VA at various stages of follow-up
# Find closest VA measurement to each of the specified snapshot times, within set tolerances
# Follow-up times in months at which the VA status is to be taken, along with acceptable tolerances
# Baseline dates must match exactly
snapshots <- data.table(months_since_index = c(0, 12, 24, 36, 48, 60, 84, 120), 
                        follow_up_month = c(0, 12, 24, 36, 48, 60, 84, 120), 
                        tolerance = c(0, 2, 2, 2, 2, 2,2,2),
                        key = "months_since_index") %>% 
  .[, `:=` (lower_lim = follow_up_month - tolerance,
            upper_lim = follow_up_month + tolerance,
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
  , `:=` (va_change_logmar = va_logmar - va_logmar[(snapshot & baseline)],
             va_change_etdrs = va_etdrs - va_etdrs[(snapshot & baseline)])
  , by=.(PatientID, EyeCode)][
  # Classify VA change into three categories
  , va_change_lines := cut(va_change_logmar, 
                               breaks = c(min(va_change_logmar), -0.2, 0.19, max(va_change_logmar)), 
                               labels = c("gained >=2 lines", "< 2 lines change", "lost >= 2 lines"), include.lowest = TRUE)]

# Note that the 'lines changed' classification breaks down for those with very low vision (e.g. counting fingers)
# va_history[PatientID=="EB9B623B-3B63-26EF-E293-64A5575590EA" & EyeCode=="L", .(EncounterDate,va_logmar, va_change_logmar, va_etdrs, va_change_etdrs, va_change_lines)]
  

    va2 %>% 
    count(va_change_logmar, va_change_etdrs, va_change_lines) %>% 
    print(n=Inf)
  
  va2 %>% filter(va_change_etdrs == 5, va_change_logmar < -0.68) 
  va2[va_change_logmar == -0.68, .(va_logmar, va_change_logmar, va_etdrs, va_change_etdrs, va_change_lines)]
  
  va2[PatientID=="3688F66F-3F74-E6FC-0276-CD9D951686A8" & EyeCode=="R", .(va_logmar, va_change_logmar, va_etdrs, va_change_etdrs, va_change_lines)]
  
  va2[PatientID=="EB2DC2F5-9BE0-4B3E-6292-9AC8D95F02BE" & EyeCode=="R", .(EncounterDate,va_logmar, va_change_logmar, va_etdrs, va_change_etdrs, va_change_lines)]
  va2[PatientID=="EB9B623B-3B63-26EF-E293-64A5575590EA" & EyeCode=="L", .(EncounterDate,va_logmar, va_change_logmar, va_etdrs, va_change_etdrs, va_change_lines)]
  
  
  va3[PatientID=="EB2DC2F5-9BE0-4B3E-6292-9AC8D95F02BE" & EyeCode=="R", .(EncounterDate,va_logmar, va_change_logmar, va_etdrs, va_change_etdrs, va_change_lines)]
  
  # Note ETDRS sometimes unreliable
  va2[PatientID=="EB2DC2F5-9BE0-4B3E-6292-9AC8D95F02BE" & EyeCode=="R" & va_change_logmar == -0.1, .(EncounterDate,va_logmar, va_change_logmar, va_etdrs, va_change_etdrs, va_change_lines)]
  
  
  
  va_history_raw[PatientID=="EB2DC2F5-9BE0-4B3E-6292-9AC8D95F02BE" & EyeCode=="R" & va_logmar == 0.1, .(EncounterDate,va_logmar, va_etdrs)]
  va_history[PatientID=="EB2DC2F5-9BE0-4B3E-6292-9AC8D95F02BE" & EyeCode=="R" & va_change_logmar == -0.1, .(EncounterDate,va_logmar, va_change_logmar, va_etdrs, va_change_etdrs, va_change_lines)]
  

  # va_history[PatientID == "F5D4AFF5-4D1B-7D14-191A-59B99B3FC8B7" & EyeCode == "R" & follow_up_month==48]
  
  va2 %>% count(va_change_etdrs, va_change_logmar, va_change_lines) %>% print(n=Inf)  

va2[va_change_etdrs == -40 & va_change_lines == "< 2 lines change",]
  
# OCT thickness history relating to the selected eyes
thickness_history <- oct_thickness_raw %>%  
  inner_join(eye %>% 
               filter(thickness_measurement) %>% 
               select(PatientID, EyeCode), by = c("PatientID", "EyeCode")) %>% 
  # Roll forward baseline thickness measurements up to 2 weeks prior to baseline
  mutate(months_since_index = if_else(baseline & months_since_index > -0.5 & months_since_index < 0, 0, months_since_index)) %>% 
  # Exclude all measurements prior to baseline
  filter(months_since_index >= 0)
  # Calculate year of treatment for each eye (since index date)
  # mutate(treatment_year = as.factor(floor(as.interval(index_date, ExamDate)/dyears()) +1 ))


# Fluid history for each eye
fluid_history <- fluid_raw %>% 
  inner_join(eye %>% 
               filter(fluid_measurement) %>% 
               select(PatientID, EyeCode), by = c("PatientID", "EyeCode")) %>% 
  # Roll forward baseline fluid measurements up to 2 weeks prior to baseline
  mutate(months_since_index = if_else(baseline & months_since_index > -0.5 & months_since_index < 0, 0, months_since_index)) %>% 
  mutate(complete_months = floor(months_since_index)) %>% 
  # Exclude all measurements prior to baseline
  filter(months_since_index >= 0) %>% 

  # Classification of AMD type by presence of SRF and IRF
  mutate(amd_type = factor(case_when(SRF == 1 & IRF == 1 ~ "SRF and IRF",
                              SRF == 0 & IRF == 1 ~ "IRF only",
                              SRF == 1 & IRF == 0 ~ "SRF only",
                              SRF == 0 & IRF == 0 ~ "No fluid"), levels = c("No fluid", "SRF only", "IRF only", "SRF and IRF")))

fluid_history %>% 
  count(amd_type, baseline)

# Find a single measurement for each eye for each month of followup




## Compare visual acuity and fluid measurements
va_fl <- eye %>% 
  inner_join(fluid_history %>% 
               filter(baseline),
             by = c("PatientID", "EyeCode"))

va_fl_correlations <- va_fl %>% 
  select(all_of(c(noa_fluid, noa_retinal, noa_grid)) & where(is.numeric)) %>% 
  map(~cor.test(.x, va_fl$va_logmar)) %>% 
  map_dfr(~broom::tidy(.), .id = "Variable") %>% 
  mutate(across(p.value, format.pval, eps = 0.001, digits = 2),
         across(c(estimate, conf.low, conf.high), .fns = ~formatC(.x, format = "f", digits = 2))) %>% 
  mutate(`Correlation with VA` = paste0(estimate, " (", conf.low, ", ", conf.high, ")")) %>% 
  transmute(Variable, estimate, `Correlation with VA`,  P = p.value)






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


# Function to produce a circular plot of the ETDRS grid with labels corresponding to values in the source data 
# Args:
# dat = data.frame with joining columns "region" and "circ" 
# label_col = bare variable name containing the labels to apply to the grid
PlotETDRSGrid <- function(dat, label_col){
  
  label_value <- enquo(label_col)
  
  etdrs_rect <- tibble(region = as_factor(c("CS", "SUP", "NAS", "INF", "TEM", "SUP", "NAS", "INF", "TEM")),
                       prefix = c("", rep("i", 4), rep("o", 4)),
                       circ = c(1, rep(3, 4), rep(6, 4)),
                       xmin = c(0, rep(0.5, 4), rep(1.5, 4)),
                       xmax = c(0.5, rep(1.5, 4), rep(3, 4)),
                       ymin = c(0, 0, 0.25, 0.5, 0.75, 0, 0.25, 0.5, 0.75),
                       ymax = c(1, 0.25, 0.5, 0.75, 1, 0.25, 0.5, 0.75,1),
                       xtext = c(0, xmax[-1] - (xmax[-1]-xmin[-1])/2),
                       ytext = c(0, ymax[-1] - (ymax[-1]-ymin[-1])/2))
  
  dat %>% 
    inner_join(etdrs_rect, by = c("region", "circ")) %>% 
  
    ggplot(aes(x = xtext, y = ytext)) +
    geom_rect(data = filter(etdrs_rect, region != "CS"), aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "white", colour = "black", size = 1.5) + 
    geom_rect(data = filter(etdrs_rect, region == "CS"), aes(xmin = xmin+0.01, xmax = xmax-0.01, ymin = ymin+0.01, ymax = ymax-0.01), fill = "white") +
    geom_text(aes(label = !!label_value)) +
    scale_y_continuous(breaks = NULL) +
    scale_x_continuous(breaks = NULL, sec.axis = sec_axis(trans = ~., breaks = NULL, name="Nasal")) +
    theme_light() +
    theme(panel.border = element_blank(), panel.grid = element_blank()) +
    labs(x= "Temporal", y = NULL)+
    coord_polar(theta = "y", start = 315*pi/180)
  
}
# 
# baseline_fluid_grid %>% 
#   filter(str_detect(abbr_metric, "IRF Vol")) %>% 
#   PlotETDRSGrid(valform) 




## Calculate treatment intervals for the selected eyes

# Just the injections relating to selected eyes
injections <- injections_clean %>% 
  inner_join(eye %>% 
               select(PatientID, EyeCode, index_date), by = c("PatientID", "EyeCode")) %>% 
  # Calculate year of treatment for each eye (since index date)
  # Starts count at year 1
  mutate(treatment_year = as.factor(floor(as.interval(index_date, EncounterDate)/dyears()) +1))

# Note that some injections for selected eyes were in other clinic types e.g. cataract or glaucoma
# This may be because an incorrect clinic type was set at the beginning of the day or because of a genuine clinical need
# These injection dates are EXCLUDED from the NOA analysis but INCLUDED in the injection interval analysis
encounters_raw %>% 
  left_join(injections, by = "EncounterID") %>% 
  count(ClinicalCategoryDesc, EncounterTypeDesc, EncounterTypeCode, InjectedDrugDesc)


treatment_intervals <- injections %>%
  group_by(PatientID, EyeCode) %>% 
  mutate(next_encounter = lead(EncounterDate),
         treatment_interval_days = interval(EncounterDate, next_encounter)/ddays(), 
         treatment_interval_weeks = interval(EncounterDate, next_encounter)/dweeks(),
         # Indicator of short treatment interval 
         treatment_interval_short = treatment_interval_days <= 40) %>% 
  filter(!is.na(next_encounter)) %>% 
  # Timeline since the index date for each person
  # First and second year since treatment start concept not relevant because do not know the date of first injection (no clinic entry date recorded)
  mutate("Treatment interval" = if_else(treatment_interval_weeks > 12, ">12 weeks", "<=12 weeks"),
         # Calculate month and year of treatment for each eye (since index date)       
         treatment_months = interval(index_date, EncounterDate)/dmonths(),
         treatment_year = as.factor(floor(as.interval(index_date, EncounterDate)/dyears()) + 1)) %>% 
  # Injection sequence (index injection = 1)
  mutate(injection_seq = row_number()) %>% 
  ungroup() %>% 
  mutate(sequence_id = abbreviate(paste(PatientID, EyeCode))) 





  




