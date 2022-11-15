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

# Fluid eyes 
fluid_eyes <- fluid_history %>% 
  count(PatientID, EyeCode, name = "measurements") %>% 
  filter(measurements > 3)


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
  mutate(cluster = fluid_kmeans[[4]]$cluster)

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

