#### Libraries ----
library(here)
library(dplyr)
library(survival)
library(survminer)
library(ggplot2)
#### Data ----

# Breast cancer data from Metabric
metabric <- readRDS(here::here("data", "metabric-analytical.rds"))


#### Functions ----


#### Pre-processing----
# Setting living as 1, dead as 0 for vital_status
metabric <-
  metabric %>%
  mutate(vital_status_level = ifelse(vital_status == "living", 0, 1))

# Factorising grade
metabric <-
  metabric %>%
  mutate(grade_level = as.factor(grade))

# Factorising cellularity
metabric <-
  metabric %>%
  mutate(cellularity_level = as.factor(cellularity))

# Factorising tumor stage
metabric <-
  metabric %>%
  mutate(tumor_stage_level = as.factor(tumor_stage))
    

#### Analysis of survival----
km <- survfit(Surv(surv_months, vital_status_level,) ~ 1, metabric)
ggsurvplot(km, title = "Survival probability over time in months", xlab = "Time in months")

#### Analysis of Logrank----
# Selecting variables for log rank
metabric_logrank <-
  metabric %>%
    select(surv_months, vital_status_level, grade_level, cellularity_level, tumor_stage_level, tumor_size)

# Logrank test for grade
grade_logrank <- 
  survdiff(Surv(surv_months, vital_status_level) ~ grade_level, data = metabric_logrank)
grade_logrank

# Logrank test for cellularity
cellularity_logrank <- 
  survdiff(Surv(surv_months, vital_status_level) ~ cellularity_level, data = metabric_logrank)
cellularity_logrank

# Logrank test for tumor stage
tumor_stage_logrank <-
  survdiff(Surv(surv_months, vital_status_level) ~ tumor_stage_level, data = metabric_logrank)
tumor_stage_logrank

# Logrank test for tumor size
tumor_size_logrank <-
  survdiff(Surv(surv_months, vital_status_level) ~ tumor_size, data = metabric_logrank)
tumor_size_logrank
  # May have to categorise them