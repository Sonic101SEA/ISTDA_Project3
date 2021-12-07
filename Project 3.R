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
# Setting died of disease as 1, rest as 0 for vital_status
metabric <- metabric %>% 
  mutate(death = ifelse(vital_status == 'died of disease', 1, 0))
  # Note: There is one missing value in vital_status

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
# Base model
km_fit_0 <- survfit(Surv(surv_months, death) ~ 1, data = metabric)
## survival prob at 60 months
summary(km_fit_0, times = c(60))
## get median
km_fit_0

## plot KM
km_fit_0_plot <- 
  ggsurvplot(km_fit_0, 
             title = "Kaplan-Meier plot for overall survivability", 
             xlab = "Time in months",
             ylab = "Overall survival probability",
             size = 0.3, 
             censor.size = 0.1,
             surv.median.line = "hv",
  )
km_fit_0_plot$plot + geom_segment(aes(x=60,xend=60, y=0, yend=0.822)) + geom_segment(aes(x=0, xend=60, y=0.822, yend=0.822))

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