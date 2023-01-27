# BIRAX - analysis of injection intervals

library(tidyverse)
library(lubridate)
library(survival)

source(rprojroot::find_rstudio_root_file("code/BIRAX-assemble-cohort.R"))

## Duration of anti-VEGF treatment ##

# All start on treatment start

eyesurv <- eye %>% 
  # Generate a single eye ID
  transmute(eye_id = paste(PatientID, EyeCode),
            Gender,
            index_age, 
            years_treated, 
            years_va, 
            years_observed, 
            PerturbedAgeatDeath,
            died = !is.na(PerturbedAgeatDeath),
            years_died = if_else(died, years_observed, as.numeric(NA)),
            # Define if treatment ending could be fully observed (i.e. not censored)
            # Define as censored if eye observed for <18 weeks after treatment end 
            observed_treatment_stop = years_observed >= years_treated + 0.35,
            # Define if end of VA followup could be fully observed
            observed_va_stop = years_observed >= years_va + 0.35
  ) 

# Convert from eye level to event level format
# Define the end of follow-up
eyesurv_mstate_raw <- tmerge(eyesurv, eyesurv, id = eye_id, tstop = years_observed) %>%
  # Split follow-up on each of the three events: treatment stop, VA monitoring stop, death
  # Treatment stop
  # If treatment stop was not fully observed (i.e. not a full 18 weeks afterwards) this gives two periods with the same state. Does not cause a modelling problem.
  tmerge(.,
    eyesurv,
    id = eye_id,
    treat_stop = event(years_treated, observed_treatment_stop)
  ) %>%
  # VA monitoring stop
  tmerge(.,
         eyesurv,
         id = eye_id,
         va_stop = event(years_va, observed_va_stop)) %>%
  # Death
  tmerge(., eyesurv, id = eye_id, death = event(years_died)) 
 
# Add row numbers for each state for each eye
eyesurv_mstate <-
  tmerge(eyesurv_mstate_raw,
         eyesurv_mstate_raw,
         id = eye_id,
         enum = cumtdc(tstart)) %>%
  
  # Mark entries after treatment stop as ended
  group_by(eye_id) %>%
  mutate(ended = if_else(is.na(dplyr::lag(treat_stop)), FALSE, dplyr::lag(cumsum(treat_stop) > 0))) %>%
  # Mark entries after VA stop as VA monitoring ended
  mutate(va_ended = if_else(is.na(dplyr::lag(va_stop)), FALSE, dplyr::lag(cumsum(va_stop) > 0))) %>%
  ungroup() %>%
  
  # Define the model event - this is the state at the end of the interval defined by times tstart and tstop
  mutate(event =
           factor(
             case_when(
               # Start with the absorbing states
               # Later definitions will not overwrite earlier definitions
               death == 1 ~ 3,
               va_stop & !treat_stop & ended ~ 2,
               va_ended & ended ~ 2,
               treat_stop ~ 1,
               ended & !va_ended ~ 1, !ended ~ 0
             )
             ,
             labels = c("Treatment", "Ended", "VA ended", "Died")
           ))

# Check distribution of states
eyesurv_mstate %>%
  count(treat_stop, va_stop, death, ended, va_ended, event)
  

injections %>% filter(PatientID == "8BE4E2CC-22B1-4C84-0C7D-ADA111C4598A", EyeCode == "R")
eye %>% filter(PatientID == "8BE4E2CC-22B1-4C84-0C7D-ADA111C4598A", EyeCode == "R")
eyesurv %>% filter(eye_id == "8BE4E2CC-22B1-4C84-0C7D-ADA111C4598A R")
eyesurv_mstate %>% filter(eye_id == "8BE4E2CC-22B1-4C84-0C7D-ADA111C4598A R")
 
 
# Examples: ended
eyesurv %>% filter(eye_id == "B0DB7775-55D0-FAD5-B328-99DB23FF3320 R")
eyesurv_mstate %>% filter(eye_id == "B0DB7775-55D0-FAD5-B328-99DB23FF3320 R")
 
# ended then died (not long enough to observe VA stop)
eyesurv %>% filter(eye_id == "C9B37F5A-CD61-D0D2-1DD8-480034E21FDB R")
eyesurv_mstate %>% filter(eye_id == "C9B37F5A-CD61-D0D2-1DD8-480034E21FDB R")
 


# Fit multistate model
msurv <-
  survfit(Surv(tstart, tstop, event) ~ 1, data = eyesurv_mstate, id = eye_id)
# Model summary
print(msurv)
# Transitions (check if plausible)
msurv$transitions
# Proportion in each state by event time (Aalen-Johanson estimate)
# Columns are in same order as rows in print(msurv) or msurv$states
msurv$pstate
mstates <- data.table(yr = msurv$time,
       `Under treatment` = msurv$pstate[,1],
       `Ended treatment` = msurv$pstate[,2],
       `VA monitoring ended` = msurv$pstate[,3],
       `Died` = msurv$pstate[,4])


# Switch to interval censoring?


### Cumulative number of anti-VEGF treatments ###

# Times of each injection in long survival format (rows for each event)
# injections_surv <- eye %>% select(PatientID, EyeCode, years_observed) %>% 
#   inner_join(injections, by = c("PatientID", "EyeCode")) %>%
#   transmute(eye_id = paste(PatientID, EyeCode),
#             years_observed,
#             tstart = years_treated) %>% 
#   group_by(eye_id) %>% 
#   mutate(tstop = lead(tstart),
#          status = 1) %>% 
#   ungroup() %>% 
#   filter(!is.na(tstop))
# Need to have end date for monitoring to show length of follow-up?
# Should be end of observation or end of VA monitoring?
# Does it make a difference?

# First tmerge is to set the overall observation period
injections_surv <- tmerge(transmute(eye, eye_id = paste(PatientID, EyeCode)),
    transmute(eye, eye_id = paste(PatientID, EyeCode), years_observed),
  id = eye_id,
  tstop = years_observed
) %>% 
  # Then add the injections as events, splitting follow-up at each event
  tmerge(.,
         transmute(injections, eye_id = paste(PatientID, EyeCode), years_treated),
         id = eye_id,
         status = event(years_treated))

injections_surv %>% head()

# Survival model for repeated event data
irep_event <- survfit(Surv(tstart, tstop, status) ~ 1, data = injections_surv, id = eye_id)

survcheck(Surv(tstart, tstop, status) ~ 1, data = injections_surv, id = eye_id)
plot(irep_event)
# Plot cumulative hazard (i.e. mean number of injections expected given the time elapsed)
plot(irep_event, cumhaz = TRUE)



# Proportions that have experienced a certain number of injections by time.

# time dependent covariate to define state
injections_states <- tmerge(transmute(eye, eye_id = paste(PatientID, EyeCode)),
              transmute(eye, eye_id = paste(PatientID, EyeCode), years_treated),
              id = eye_id,
              tstart = -0.001,
              tstop = years_treated) %>% 
  tmerge(.,
         injections %>% 
           transmute(eye_id = paste(PatientID, EyeCode), years_treated) %>% 
           group_by(eye_id) %>% 
           # Initial state is 1 (all have an injection at years_treated=0)
           mutate(injstate = row_number()),
         id = eye_id,
         status = event(years_treated, injstate)) %>% 
  # State must be a factor for multi-state modelling
  mutate(istate = as.factor(status))

# Loses the time zero states.

injections %>% head()
injections_states %>% filter(eye_id == "000EFFF4-3FD0-760C-BB88-84BA43222F23 L")

injections_states %>% count(istate)

injections_states %>% 
  count(eye_id) %>% 
  filter(n==3) %>% 
  inner_join(injections_states)

survcheck(Surv(tstart, tstop, istate) ~ 1, data = injections_states, id = eye_id)
# Slow
isurv <- survfit(Surv(tstart, tstop, istate) ~ 1, data = injections_states, id = eye_id)

isurv$transitions

plot(isurv)
# Plot cumulative hazard (i.e. mean number of injections expected given the time elapsed)
plot(isurv, cumhaz = TRUE)


# Numbers of injections received (categorical) by time
istates <- data.table(yr = isurv$time,
                      isurv$pstate)
iarea <- istates %>% 
  pivot_longer(-yr, names_to = "State", values_to = "Proportion in state") %>% 
  mutate(State_num = as.numeric(substr(State, 2, nchar(State))),
           State = fct_rev(as_factor(State_num))) %>% 
  mutate(n_injection = fct_rev(cut(State_num, breaks = c(1, 2, 3, 4, 5, 10, 20, 50, max(.$State_num)), include.lowest = TRUE)))




