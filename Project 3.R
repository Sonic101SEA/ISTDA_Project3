#### Libraries ----
library(dplyr)
library(tidyverse)
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

## Descriptive analysis of variables
# For receptors
long_table_receptors <-
  metabric %>%
  pivot_longer(cols = c(er_status, her2_status, pr_status),
               names_to = "variable",
               values_to = "value") %>%
  group_by(variable, value) %>%
  summarise(n = n()) %>%
  mutate(value = factor(
    value,
    levels = c("negative", "positive")))

barplot_receptors <- long_table_receptors %>%
  ggplot(aes(variable, n, fill = value)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Distribution of positive and negative receptor status", 
       x = "Receptor Status", y = "Count", fill = "Status")

# Without using group_by and summarise
long_table_receptors2 <-
  metabric %>%
  pivot_longer(cols = c(er_status, her2_status, pr_status),
               names_to = "variable",
               values_to = "value") %>%
  mutate(value = factor(
    value,
    levels = c("negative", "positive")))

barplot_receptors2 <- long_table_receptors2 %>%
  ggplot(aes(variable, fill = value)) +
  geom_bar(stat = "count") +
  labs(title = "Distribution of positive and negative receptor status", 
       x = "Receptor Status", y = "Count", fill = "Status")

# Individual bar plots for receptors
barplot_er <- metabric %>%
  ggplot(aes(x = er_status)) +
  geom_bar() +
  geom_text(stat = "count", aes(label=..count..), vjust=-1)

barplot_her2 <- metabric %>%
  ggplot(aes(x = her2_status)) +
  geom_bar()

barplot_pr <- metabric %>%
  ggplot(aes(x = pr_status)) +
  geom_bar()

# Individual barplot for treatment
long_table_treatment <-
  metabric %>%
  pivot_longer(cols = c(chemotherapy, hormone_therapy, radio_therapy),
               names_to = "variable",
               values_to = "value") %>%
  group_by(variable, value) %>%
  summarise(n = n()) %>%
  mutate(value = factor(
    value,
    levels = c("no", "yes")))

barplot_treatment <- long_table_treatment %>%
  ggplot(aes(variable, n, fill = value)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Distribution of treatment status", 
       x = "Treatment Status", y = "Count", fill = "Status")

# Barplot of vital_status
long_table_vital <-
  metabric %>%
  pivot_longer(cols = c(vital_status),
               names_to = "variable",
               values_to = "value") %>%
  group_by(variable, value) %>%
  summarise(n = n()) %>%
  mutate(value = factor(
    value))

barplot_vital <- long_table_vital %>%
  ggplot(aes(value, n, fill = n, label = n)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Distribution of Vital status", 
       x = "Vital Status", y = "Count", fill = "Status") +
  geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  theme(legend.position = "none")

# Boxplot for age_diagnosis
metabric %>%
  ggplot(aes(x = age_diagnosis, y = "Age at diagnosis")) + 
  geom_boxplot(width = 0.2) +
  labs(title = "Distribution of Age at diagnosis",
       x = "Age in months") +
  theme(axis.title.x = element_blank()) +
  coord_flip()

# Boxplot for surv_months
metabric %>%
  ggplot(aes(x = surv_months, y = "Survival")) + 
  geom_boxplot(width = 0.2) +
  labs(title = "Distribution of survival months",
       x = "Survival in months") +
  theme(axis.title.x = element_blank()) +
  coord_flip()

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
km_fit_0_plot$plot + 
  geom_segment(aes(x=60,xend=60, y=0, yend=0.822), colour = "grey60", linetype = "dotted") + 
  geom_segment(aes(x=0, xend=60, y=0.822, yend=0.822), colour = "grey60", linetype = "dotted") +
  annotate("text", x = 200, y = 0.822, label = "5-year survival probability = 0.822") +
  annotate("text", x = 320, y = 0.60, label = "Median survival = 283 months")

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