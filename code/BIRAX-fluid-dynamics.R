# BIRAX - analysis of retinal fluid dynamics
library(tidyverse)
library(lubridate)
library(survival)
library(parallelDist)

source(rprojroot::find_rstudio_root_file("code/BIRAX-assemble-cohort.R"))

# Probability of fluid presence by month of follow-up

# Only include months with a certain number of measurements to ensure stability of estimates
trunc_n <- 100
trunc_month <- fluid_history %>% 
  count(treatment_months) %>% 
  filter(n > trunc_n) %>% 
  slice_min(n) %>% 
  pull(treatment_months)
  

# Calculate quantiles of fluid distribution by treatment month
fluid_qtl <- fluid_history %>% 
  group_by(treatment_months) %>% 
  nest() %>% 
  mutate(irf_qtl = map(data, ~enframe(quantile(.$IRFVolumeNl, probs = c(0.5, 0.6, 0.7, 0.8, 0.9), na.rm = TRUE), name = "Quantile", value = "IRF")),
         srf_qtl = map(data, ~enframe(quantile(.$SRFVolumeNl, probs = c(0.5, 0.6, 0.7, 0.8, 0.9), na.rm = TRUE), name = NULL, value = "SRF"))) %>% 
  unnest(cols = c(irf_qtl, srf_qtl)) %>% 
  pivot_longer(cols = c(IRF, SRF), names_to = "fluid_type", values_to = "volume")


# Dynamic time warping

# Extract the multivariate fluid sequences for each eye
fluid_seq_raw <- fluid_history %>% 
  inner_join(fluid_eyes, by = c("PatientID", "EyeCode")) %>% 
  group_by(PatientID, EyeCode) %>% 
  summarise(IRFVolumeNl, SRFVolumeNl, ti = row_number(), .groups = "keep") %>% 
  nest() %>% 
  mutate(fseq = map(data, ~t(select(., IRFVolumeNl, SRFVolumeNl))))

# Align the sequences and compute the distances between the aligned sequences
fluid_dist <- parDist(fluid_seq_raw$fseq, method = "dtw", open.end = TRUE)

# Cluster by distance
# Run KMeans for different numbers of clusters
# 2 mins
fluid_kmeans <- vector("list", length = 10)
for(i in 1:length(fluid_kmeans)){
  fluid_kmeans[[i]] <- kmeans(fluid_dist, centers = i)   
}

# Add the cluster ID for each eye
fluid_seq <- fluid_seq_raw %>% 
  ungroup() %>% 
  # k selected using scree plot
  mutate(cluster = as.factor(fluid_kmeans[[4]]$cluster))

# Subset fluid eye and VA tables for just those eyes in the cluster analysis
eye_cluster <- eye %>% 
  inner_join(fluid_seq %>% 
               select(PatientID, EyeCode, cluster), 
             by = c("PatientID", "EyeCode"))

fluid_history_cluster <- fluid_history %>% 
  inner_join(eye_cluster %>% 
               select(PatientID, EyeCode, cluster),
             by = c("PatientID", "EyeCode"))

va_history_cluster <- va_history %>% 
  inner_join(eye_cluster %>% 
               select(PatientID, EyeCode, cluster),
             by = c("PatientID", "EyeCode"))

thickness_history_cluster <- thickness_history %>% 
  inner_join(eye_cluster %>% 
               select(PatientID, EyeCode, cluster),
             by = c("PatientID", "EyeCode"))


# Calculate quantiles for each of the clusters
fluid_qtl_cluster <- fluid_seq %>% 
  select(PatientID, EyeCode, cluster) %>% 
  inner_join(fluid_seq %>% 
               count(cluster), by = "cluster") %>% 
  inner_join(fluid_history, by = c("PatientID", "EyeCode")) %>% 
  group_by(treatment_months, cluster, n) %>% 
  nest() %>% 
  mutate(irf_qtl = map(data, ~enframe(quantile(.$IRFVolumeNl, probs = c(0.5, 0.6, 0.7, 0.8, 0.9), na.rm = TRUE), name = "Quantile", value = "IRF")),
         srf_qtl = map(data, ~enframe(quantile(.$SRFVolumeNl, probs = c(0.5, 0.6, 0.7, 0.8, 0.9), na.rm = TRUE), name = NULL, value = "SRF"))) %>% 
  unnest(cols = c(irf_qtl, srf_qtl)) %>% 
  pivot_longer(cols = c(IRF, SRF), names_to = "fluid_type", values_to = "volume") %>% 
  ungroup()


# GLMs for differences in descriptive statistics by cohort
summary(lm(years_observed ~ cluster, data = eye_cluster))
summary(lm(years_treated ~ cluster, data = eye_cluster))
summary(lm(total_injections ~ cluster, data = eye_cluster))


# Longitudinal outcomes in terms of visual acuity

# How to measure? 
# VA below a certain threshold?
# Event at earliest drop below threshold or lose two lines?

# Calculate quantiles for each of the clusters
va_qtl_cluster <- va_history_cluster %>% 
  inner_join(fluid_seq %>% 
               count(cluster), by = "cluster") %>%  
  group_by(treatment_months, cluster, n) %>% 
  nest() %>% 
  mutate(logmar_qtl = map(data, ~enframe(quantile(.$va_logmar, probs = c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm = TRUE), name = "Quantile", value = "va_logmar"))) %>% 
  unnest(cols = logmar_qtl) %>% 
  ungroup()




# Predicting cluster membership using baseline characteristics


# Multiclass classification problem

m1 <- nnet::multinom(cluster ~ Gender + va_logmar + IRFVolumeNl, data = eye_cluster %>% inner_join(filter(fluid_history_cluster, baseline)))
summary(m1)

broom::tidy(m1) %>% 
  mutate(p.value = format.pval(p.value, eps = 0.001))

table(predict(m1))

# ROC performance of IRF volume alone
fluid_history_cluster %>% 
filter(baseline) %>% 
  mutate(cluster3 = cluster == 3) %>% 
  pROC::roc(data = ., response = cluster3, predictor = IRFVolumeNl) %>% 
  plot()


# Imputation of missing variables

noa_fluid
noa_general

# All volume measurements except total volumes, as these are linear combinations of IRF and SRF
noa_volumes <- noa %>% 
  select(matches("volume|_vol"), -matches("TRF")) %>% 
  names()

# noa %>% 
#   filter(IRFVolumeNl + SRFVolumeNl != TRFVolumeNl) %>% 
#   mutate(IRFVolumeNl + SRFVolumeNl) %>% 
#   select(matches("VolumeNl"))

# Quantify missingness for volumes
fluid_history_cluster %>% 
  select(all_of(noa_volumes)) %>% 
  summarise(across(.fns = ~sum(is.na(.))/length(.)*100)) %>% 
  pivot_longer(cols = everything()) %>% 
  arrange(desc(value)) %>% print(n=Inf)

# Just select volumes with a negligible amount of imputation
noa_volumes_selected <- noa %>% 
  select(all_of(noa_volumes), -matches("Superior6|Inferior6")) %>% 
  names()

noa_volumes_selected <- noa_volumes

# Impute missing volumes with zero and scale
fluid_history_imp <- fluid_history_cluster %>% 
  mutate(across(all_of(noa_volumes), ~scale(if_else(is.na(.), 0, .)))) %>% 
  arrange(PatientID, EyeCode, months_since_index) %>% 
  group_by(PatientID, EyeCode) %>% 
  mutate(oct_series = row_number()) %>% 
  ungroup() %>% 
  as.data.table()


fluid_baseline_imp <- fluid_history_imp %>% 
  select(cluster, PatientID, EyeCode, oct_series, all_of(noa_volumes_selected)) %>% 
  mutate(cluster = as.numeric(cluster)) %>% 
  filter(oct_series == 1)


# Just those with a snapshot at 6 months
fluid_6months_imp <- fluid_history_imp %>% 
  filter(snapshot, follow_up_month <= 6) %>% 
  select(cluster, PatientID, EyeCode, follow_up_month, all_of(noa_volumes_selected)) %>% 
  mutate(cluster = as.numeric(cluster)) %>% 
  pivot_wider(id_cols = c(cluster, PatientID, EyeCode), 
              names_from = follow_up_month,
              values_from = all_of(noa_volumes_selected)) %>% 
  na.omit()

# Just those with a snapshot at 6 months and 12 months
fluid_6_12_imp <- fluid_history_imp %>% 
  filter(snapshot, follow_up_month <= 12) %>% 
  select(cluster, PatientID, EyeCode, follow_up_month, all_of(noa_volumes_selected)) %>% 
  mutate(cluster = as.numeric(cluster)) %>% 
  pivot_wider(id_cols = c(cluster, PatientID, EyeCode), 
              names_from = follow_up_month,
              values_from = all_of(noa_volumes_selected)) %>% 
  na.omit()


fluid_6months_imp %>% complete.cases()

fluid_history_imp %>% count(follow_up_month, snapshot)

# Drop NOA measurements ruled ineligible for analysis
fluid_history_cluster %>% 
  count(`Analysis eligibility`)


# Calculate change variables?

